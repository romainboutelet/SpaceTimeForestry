rm(list=ls())
library(tidyverse)
library(sf)
library(FIESTA)
library(spdep)
library(dplyr)
library(pals)

# Load in CONUS shapefile from FIESTA package
shp <- stunitco[!(stunitco$STATENM %in% c("Alaska", "Hawaii", "Puerto Rico", "Commonwealth of the Nort*",
                                        "Guam", "United States Virgin Isl*",
                                        "American Samoa")), c("STATENM", "COUNTYNM", "COUNTYFIPS")]

colnames(shp) <- c("state", "county", "fips", "geometry")
# shp$fips <- as.numeric(shp$fips)
# Add TCC data (row order is already matching shp)
tcc <- read.csv("../mean_tcc_fia_county_level_2008_2021.csv")
shp <- cbind(shp, tcc[, -1])

eco <- st_read("../na_cec_eco_l2/NA_CEC_Eco_Level2.shp") %>% st_transform(crs = st_crs(shp)) %>% filter(NA_L2NAME != "WATER")
shp <- cbind(shp, "eco" = str_to_title(eco[st_nearest_feature(st_centroid(shp), eco),]$NA_L2NAME))

# Load in FIA data
dat <- readRDS("../carbon_density_on_all_plots_CONUS_measyear2000_to_present_w_CWD.Rds")

years <- 2008:2021
dat <- dat %>% filter(measyear %in% years)

dat$fips <- paste0(sprintf("%02d", dat$statecd), sprintf("%03d", dat$countycd))
dat$id <- paste0(sprintf("%02d", dat$statecd), sprintf("%01d", dat$unitcd),
                 sprintf("%03d", dat$countycd), sprintf("%05d", dat$plot))

dat <- dat %>% filter(fips %in% st_drop_geometry(shp$fips))

dat <- dat %>% filter(intensity == 1)

fips_to_eco <- st_drop_geometry(shp)[,c("fips", "eco")]

dat <- dat %>% mutate(live_trees_tons_c_per_acre = if_else(is.na(ag_live_tree_c_tons_per_acre), 0, ag_live_tree_c_tons_per_acre),
                      dead_trees_tons_c_per_acre = if_else(is.na(ag_dead_tree_c_tons_per_acre), 0, ag_dead_tree_c_tons_per_acre),
                      cwd_tons_c_per_acre = cond_CWD_c_tons_per_acre,
                      live_trees_c_Mg_ha = 2.24170023*live_trees_tons_c_per_acre,
                      dead_trees_c_Mg_ha = 2.24170023*dead_trees_tons_c_per_acre,
                      cwd_c_Mg_ha = 2.24170023*cwd_tons_c_per_acre)
                      
# Correct two county fips that are incorrect
dat[dat$fips == "12025",]$fips <- "12086"
dat[dat$fips == "46113",]$fips <- "46102"

dat <- dat %>% left_join(fips_to_eco, by = "fips")

paste0(shp[which(!(shp$fips %in% dat$fips)), c("state", "county")]$county, ", ",shp[which(!(shp$fips %in% dat$fips)), c("state", "county")]$state)

ggplot() +
  geom_sf(data = shp[shp$state == "Virginia",], aes(fill = n_distinct), col = "white", lwd = 0.05) +
  scale_fill_gradientn(colors = ocean.haline(n = 100)) + 
  geom_sf(data = shp[which(!(shp$fips %in% dat$fips) & shp$state != "Colorado"),], fill = "red", col = "white") + 
  theme_void()

ggplot() +
  geom_sf(data = shp[shp$state == "Illinois",], aes(fill = n_distinct), col = "white", lwd = 0.2) +
  scale_fill_gradientn(colors = ocean.haline(n = 100)) + 
  theme_void()

plots <- dat %>% dplyr::select(id, forest_prop, fips, eco, measyear,
                       "live_Mg_ha" = live_trees_c_Mg_ha,
                       "dead_Mg_ha" = dead_trees_c_Mg_ha,
                       "cwd_Mg_ha" = cwd_c_Mg_ha)

conus <- shp

plots[plots$forest_prop == 0, "cwd_Mg_ha"] <- 0

duplicates <- which(duplicated(plots[,c("id", "measyear")]) | duplicated(plots[,c("id", "measyear")], fromLast = TRUE))

plots <- plots[-duplicates[seq(from = 1, to = length(duplicates), by = 2)],]

n_jt <- plots %>% dplyr::group_by(fips, measyear) %>% summarize(n_jt = n())

no_plots <- which(!(shp$fips %in% dat$fips))

for (t in years) {
  n_t <- n_jt[n_jt$measyear == t, c("fips", "n_jt")]
  colnames(n_t) <- c("fips", paste0("n_", t))
  conus <- left_join(conus, n_t, by = "fips") 
}

for (t in years) {
  conus[which(is.na(conus[,paste0("n_", t)])), paste0("n_", t)] <- 0
}

M <- 3
P <- 2 # predictors including intercept
Q <- 1 # Space-varying predictors no intercept
J <- nrow(conus)
T <- length(years)
n_jt <- st_drop_geometry(conus[,c("fips", paste0("n_", years))])
colnames(n_jt) <- c("fips", years)
regions <- unique(st_drop_geometry(conus$eco))
R <- length(regions)
j_r <- match(st_drop_geometry(conus$eco), regions)

y <- list()

for (t in 1:T) {
  y_t <- list()
  for (j in 1:J) {
    fip <- st_drop_geometry(conus$fips)[j]
    y_j <- matrix(data = NA, nrow = M, ncol = n_jt[n_jt$fips == fip, paste0(years[t])],
                  dimnames = list(c("live", "dead", "cwd"), 
                                  unname(unlist(c(plots[plots$fips == fip & plots$measyear == years[t], "id"])))))
    y_t[[fip]] <- y_j
  }
  y[[as.character(years[t])]] <- y_t
}

for (i in 1:nrow(plots)) {
  year <- as.character(plots[i, "measyear"])
  fip <- as.character(plots[i, "fips"])
  id <- as.character(plots[i, "id"])
  y[[year]][[fip]][,id] <- as.numeric(plots[i, c(6,7,8)])
}

X <- array(NA, dim = c(M, M*P, J, T))
X_tilde <- array(NA, dim = c(M, M, J, T))
for (t in 1:T) {
  for (j in 1:J) {
    X[,,j,t] <- diag(1, nrow = M) %x% t(unlist(c(1, unname(st_drop_geometry(conus[j, paste0("tc_", years[t])])))))
    X_tilde[,,j,t] <- diag(unname(st_drop_geometry(conus[j, paste0("tc_", years[t])])), nrow = M)
  }
}

links <- poly2nb(conus)
W <- matrix(0, nrow = J, ncol = J)
for (j in 1:J) {
  W[j, links[[j]]] <- 1
}

W[which(shp$county == "San Juan" & shp$state == "Washington"),
  which(shp$county %in% c("Whatcom", "Skagit") & shp$state == "Washington")] <- 1

W[which(shp$county %in% c("Whatcom", "Skagit") & shp$state == "Washington"),
  which(shp$county == "San Juan" & shp$state == "Washington")] <- 1

W[which(shp$county == "Island" & shp$state == "Washington"),
  which(shp$county %in% c("San Juan", "Skagit", "Snohomish", "Kitsap", "Jefferson") & shp$state == "Washington")] <- 1

W[which(shp$county %in% c("San Juan", "Skagit", "Snohomish", "Kitsap", "Jefferson") & shp$state == "Washington"),
  which(shp$county == "Island" & shp$state == "Washington")] <- 1

W[which(shp$county == "Nantucket" & shp$state == "Massachusetts"),
  which(shp$county %in% c("Dukes", "Barnstable") & shp$state == "Massachusetts")] <- 1

W[which(shp$county %in% c("Dukes", "Barnstable") & shp$state == "Massachusetts"),
  which(shp$county == "Nantucket" & shp$state == "Massachusetts")] <- 1

D <- diag(rowSums(W))

data_list <- list("shp" = conus,
                  "plots" = plots,
                  "y" = y,
                  "X" = X,
                  "X_tilde" = X_tilde,
                  "W" = W,
                  "D" = D,
                  "years" = years,
                  "M" = M,
                  "P" = P, # predictors including intercept
                  "Q" = Q, # Space-varying predictors no intercept
                  "J" = J,
                  "T" = T,
                  "n_jt" = n_jt,
                  "regions" = regions,
                  "R" = R,
                  "j_r" = j_r
                  )

saveRDS(data_list, "./output/plot_level_carbon_pools_conus.rds")
