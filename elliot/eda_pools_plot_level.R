library(sf)
library(ggplot2)
library(dplyr)
library(pals)
library(patchwork)

test <- readRDS("./output/plot_level_carbon_pools.rds")

plots <- test[[1]]
shp <- test[[2]]

ecoregions <- unique(plots$eco)
for (e in ecoregions) {
  plots_e <- plots[plots$eco == e,]
  shp_e <- st_drop_geometry(shp[shp$eco == e,])
  plots_e$tcc <- NA
  for(i in 1:nrow(plots_e)) {
    plots_e[i, "tcc"] <- shp_e[which(shp_e$fips == as.character(plots_e[i, "fips"])), paste0("tc_", plots_e[i, "measyear"])]
  }
  
  lm_live_e <- lm(plots_e$live_Mg_ha ~ plots_e$tcc)
  lm_dead_e <- lm(plots_e$dead_Mg_ha ~ plots_e$tcc)
  lm_cwd_e <- lm(plots_e$cwd_Mg_ha ~ plots_e$tcc)
  
  live_not_na <- which(!is.na(plots_e$live_Mg_ha))
  dead_not_na <- which(!is.na(plots_e$dead_Mg_ha))
  cwd_not_na <- which(!is.na(plots_e$cwd_Mg_ha))
  
  resids <- matrix(NA, nrow = nrow(plots_e), ncol = 3)
  colnames(resids) <- c("live", "dead", "cwd")
  resids[live_not_na, "live"] <- unname(residuals(lm_live_e))
  resids[dead_not_na, "dead"] <- unname(residuals(lm_dead_e))
  resids[cwd_not_na, "cwd"] <- unname(residuals(lm_cwd_e))
  
  print(e)
  print(cor(resids, use = "pairwise.complete.obs"))
}

means <- plots %>% filter(measyear == 2013) %>% group_by(fips) %>% summarize("m_live" = mean(live_Mg_ha, na.rm = TRUE),
                                                                                  "m_dead" = mean(dead_Mg_ha, na.rm = TRUE),
                                                                                  "m_cwd" = mean(cwd_Mg_ha, na.rm = TRUE))

missing_fips <- setdiff(shp$fips, means$fips)

means <- rbind(means, data.frame("fips" = missing_fips,
                                 "m_live" = rep(NA, times = length(missing_fips)),
                                 "m_dead" = rep(NA, times = length(missing_fips)),
                                 "m_cwd" = rep(NA, times = length(missing_fips))))


means$fips <- factor(means$fips, levels = shp$fips)

means <- means %>% arrange(fips)

layout <- "
AABB
AABB
#CC#
#CC#
"

ggplot() + 
  geom_sf(data = shp, aes(fill = means$m_live), col = "grey", lwd = 0.01) +
  scale_fill_gradientn(colors = rev(ocean.algae(n = 100))) +
  labs(fill = "Live C\nMg/ha") +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  ggplot() + 
  geom_sf(data = shp, aes(fill = means$m_dead), col = "grey", lwd = 0.01) +
  scale_fill_gradientn(colors = rev(ocean.turbid(n = 100))) +
  labs(fill = "Dead C\nMg/ha") +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  ggplot() +
  geom_sf(data = shp, aes(fill = means$m_cwd), col = "grey", lwd = 0.01) +
  scale_fill_gradientn(colors = ocean.amp(n = 100)) +
  labs(fill = "CWD C\nMg/ha") +
  theme_void() +
  theme(text = element_text(family = "serif")) +
  plot_layout(design = layout) +
  plot_annotation(
    title = "2013",
    theme = theme(
      plot.title = element_text(size = 20, family = "serif", hjust = 0.5)
    )
  )
