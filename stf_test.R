rm(list=ls())
library(tidyverse)
library(sf)
library(leaflet)
library(FIESTA)
library(spdep)
library(dplyr)
library(pals)


# Add data (row order is already matching shp)
dat <- read.csv("Middle Rockies(in).csv")

# Load in CONUS shapefile from FIESTA package
dat$fips <- paste0(sprintf("%02d", dat$statecd), sprintf("%03d", dat$countycd))

dat$id <- paste0(sprintf("%02d", dat$statecd),
                 sprintf("%01d", dat$unitcd),
                 sprintf("%03d", dat$countycd),
                 sprintf("%05d", dat$plot))

shp <- stunitco[stunitco$COUNTYFIPS %in% dat$fips,]

sf_dat <- st_as_sf(dat, coords = c("lon_public", "lat_public"), crs = 5070)


my_sf <- ggplot() + geom_sf(data = subset(sf_dat, fips == 16043),
                            aes(color = measyear)) +
  stat_sf_coordinates(data = subset(sf_dat, drybio_ag_live_5plus == 0 & 
                                      fips == 16043),
                      aes(), col = "red") +
  scale_color_viridis_c()

sf_data_wgs84 <- st_transform(sf_dat, crs = 4326)

sf_data_wgs84$colors <- "blue"
sf_data_wgs84$colors[sf_data_wgs84$drybio_ag_live_5plus == 0] <- "red"

leaflet(subset(sf_data_wgs84, id != 0)) %>%
  # Add an Esri basemap (e.g., Esri.WorldImagery)
  addProviderTiles(providers$Esri.WorldImagery) %>%
  # Add your sf object to the map (default style for sf objects)
  addCircleMarkers(
    radius = ~(drybio_ag_live_5plus)/2,
    color = ~colors
    )

new_count <- dat %>% 
  group_by(id) %>% 
  summarise("n" = n())

newlist <- list()
count <- 1
for (i in subset(new_count, id %in% new_count[new_count$n > 2,]$id)$id){
  mydat <- subset(dat, id == i)
  for (j in 1:(nrow(mydat)-1)){
    newlist[[count]] <- data.frame(id = i, measyear = mydat$measyear[j],
                                   diff = mydat$drybio_ag_live_5plus[j+1]-
                                     mydat$drybio_ag_live_5plus[j],
                                   range = max(mydat$drybio_ag_live_5plus)-
                                     min(mydat$drybio_ag_live_5plus))
    count <- count + 1
  }
}

worst_dat <- do.call(rbind, newlist)

worst <- which.min(worst_dat$diff)



my_id <- subset(worst_dat, id %in% new_count[new_count$n > 2,]$id)[worst,]$id


for (i in subset(dat, id %in% new_count[new_count$n > 2,]$id)$id){
  x <- subset(dat, id == i)$measyear
  y <- subset(dat, id == i)$drybio_ag_live_5plus
  plot(x[order(x)], y[order(x)],
       type = "l", ylim = range(0, subset(dat, id == i)$drybio_ag_live_5plus))
  readline(prompt = "Pause. Press <Enter> to continue...")
}


sf_blackhills <- st_crop(sf_data_wgs84, xmin = -106, xmax = -90, ymin = 0,
                         ymax = 90)

leaflet(subset(sf_blackhills, id != 0)) %>%
  # Add an Esri basemap (e.g., Esri.WorldImagery)
  addProviderTiles(providers$Esri.WorldImagery) %>%
  # Add your sf object to the map (default style for sf objects)
  addCircleMarkers(
    radius = ~(drybio_ag_live_5plus)/2,
    color = ~colors
  )


new_count_bh <- sf_blackhills %>% 
  group_by(id) %>% 
  summarise("n" = n())

newlist_bh <- list()
count <- 1
for (i in subset(sf_blackhills, id %in% new_count_bh[new_count_bh$n > 2,]$id)$id){
  mydat <- subset(sf_blackhills, id == i)
  for (j in nrow(mydat)-1){
    newlist_bh[[count]] <- data.frame(id = i, measyear = mydat$measyear[j],
                                   diff = mydat$drybio_ag_live_5plus[j+1]-
                                     mydat$drybio_ag_live_5plus[j],
                                   range = max(mydat$drybio_ag_live_5plus)-
                                     min(mydat$drybio_ag_live_5plus))
    count <- count + 1
  }
}

worst_sf_bh <- do.call(rbind, newlist_bh)

worst <- which.max(worst_dat$range)



my_id <- subset(dat, id %in% new_count[new_count$n > 2,]$id)[worst,]$id


for (i in subset(dat, id %in% new_count[new_count$n > 2,]$id)$id){
  x <- subset(dat, id == i)$measyear
  y <- subset(dat, id == i)$drybio_ag_live_5plus
  plot(x[order(x)], y[order(x)],
       type = "l", ylim = range(0, subset(dat, id == i)$drybio_ag_live_5plus))
  readline(prompt = "Pause. Press <Enter> to continue...")
}



