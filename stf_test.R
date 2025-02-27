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

new_dat <- dat %>% 
  complete(measyear, id) %>%
  arrange(id, measyear) %>%
  group_by(id) %>%
  fill(plt_cn, statecd, unitcd, countycd, plot, drybio_ag_live_5plus, 
       lat_public, lon_public, puid, measmon, measday, time, fips) %>%
  mutate(deforested = NA, reforested = NA, year_def = NA, year_ref = NA) %>%
  fill(plt_cn, statecd, unitcd, countycd, plot, lat_public, lon_public, puid,
       measmon, measday, time, fips, .direction = "up")

for (measyear in 2002:2021){
  new_dat[new_dat$measyear == measyear,]$reforested <- 
    (new_dat[new_dat$measyear == measyear-1,]$drybio_ag_live_5plus == 0) &
    (new_dat[new_dat$measyear == measyear,]$drybio_ag_live_5plus > 0)
  new_dat[new_dat$measyear == measyear,]$deforested <- 
    (new_dat[new_dat$measyear == measyear-1,]$drybio_ag_live_5plus > 0) &
    (new_dat[new_dat$measyear == measyear,]$drybio_ag_live_5plus == 0)
}



for (measyear in 2021:2002){
  mask <- (new_dat[new_dat$measyear == measyear-1,]$drybio_ag_live_5plus > 0) &
    (new_dat[new_dat$measyear == measyear,]$drybio_ag_live_5plus == 0)
  new_dat[new_dat$measyear == measyear,]$year_def[mask] <- measyear
  mask <- (new_dat[new_dat$measyear == measyear-1,]$drybio_ag_live_5plus == 0) &
    (new_dat[new_dat$measyear == measyear,]$drybio_ag_live_5plus > 0)
  new_dat[new_dat$measyear == measyear,]$year_ref[mask] <- measyear
}


new_dat <- new_dat %>% mutate(color = "grey")
new_dat[which(new_dat$deforested), "color"] <- "red"
new_dat[which(new_dat$reforested), "color"] <- "green"

sf_new_dat <- st_as_sf(new_dat, coords = c("lon_public", "lat_public"), crs = 5070)

sf_new_data_wgs84 <- st_transform(sf_new_dat, crs = 4326)

pal <- colorNumeric(palette = rainbow(15), domain = 2006:2021)

leaflet(subset(sf_new_data_wgs84, deforested)) %>%
  # Add an Esri basemap (e.g., Esri.WorldImagery)
  addProviderTiles(providers$Esri.WorldImagery) %>%
  # Add your sf object to the map (default style for sf objects)
  addCircleMarkers(
    radius = ~5,
    color = ~pal(year_def)
  )

2006 "#FF0000"
2007 "#FF6200"
2008 "#FFC000"
2009 "#D8F500"
2010 "#87FF00"
2011 "#3AFF00"
2012 "#0BFF4C"
2013 "#20FF9E"
2014 "#21E7E4"
2015 "#28A3FF"
2016 "#134FFF"
2017 "#3300FF"
2018 "#7F00FF"
2019 "#D400F8"
2020 "#FF00C5"
2021 "#FF0066"

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
sf_data_wgs84$colors[sf_data_wgs84$id == "46208120233"] <- "green"


leaflet(subset(sf_data_wgs84, id %in% subset(new_count, n == 3)$id)) %>%
  # Add an Esri basemap (e.g., Esri.WorldImagery)
  addProviderTiles(providers$Esri.WorldImagery) %>%
  # Add your sf object to the map (default style for sf objects)
  addCircleMarkers(
    radius = ~1,
    color = ~colors
    )

new_count <- dat %>% 
  group_by(id) %>% 
  summarise("n" = n()) %>%
  mutate("zero" = NA)

for (i in 1:nrow(new_count)){
  new_count$zero[i] <- 
    (any(subset(dat, id == new_count[i,]$id)$drybio_ag_live_5plus == 0) & 
    !all(subset(dat, id == new_count[i,]$id)$drybio_ag_live_5plus == 0))
}

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

## Plotting plots with at least 2 values and that have at least one 0 but not only 0's

for (i in subset(dat, id %in% new_count[(new_count$n > 2) & new_count$zero,]$id)$id){
  x <- subset(dat, id == i)$measyear
  y <- subset(dat, id == i)$drybio_ag_live_5plus
  plot(x[order(x)], y[order(x)],
       type = "l", ylim = range(0, subset(dat, id == i)$drybio_ag_live_5plus))
  print(i)
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



