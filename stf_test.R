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

dat_count <- dat %>%
  group_by(id) %>%
  summarise("n_count" = n())

dat <- left_join(dat, dat_count, by = id)

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

leaflet(subset(sf_data_wgs84, id == "16308186024")) %>%
  # Add an Esri basemap (e.g., Esri.WorldImagery)
  addProviderTiles(providers$Esri.WorldImagery) %>%
  # Add your sf object to the map (default style for sf objects)
  addCircleMarkers(
    radius = ~(drybio_ag_live_5plus)/5,
    color = ~colors
    )

subset(dat, n_count > 2)

