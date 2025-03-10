library(terra)
library(sf)
library(stars)
library(FIESTA)
library(dplyr)
library(ggplot2)
library(exactextractr)

# Create shape file

str_folder <- '/home/romain/Documents/Research/fia_tcc_CONUS_all'

# Load in CONUS shapefile from FIESTA package
shp <- stunitco[(stunitco$STATENM == "Delaware"), c("STATENM", "COUNTYNM", "COUNTYFIPS")] %>%
  st_union()

plot(shp)

st_write(shp, paste0(str_folder, "/delaware.shp"), driver = "ESRI Shapefile")

## Need to do this in terminal:
#gdalwarp -cutline delaware.shp -crop_to_cutline -dstalpha science_tcc_CONUS_2019_v2021-4/science_tcc_conus_2019_v2021-4.tif delaware_tcc.tif

# Exctracting TCC data

#str_name <- "/science_tcc_CONUS_2019_v2021-4/science_tcc_conus_2019_v2021-4.tif"
str_name <- "/delaware_tcc.tif"

myraster <- rast(paste0(str_folder,str_name))

tcc_del <- myraster[myraster$delaware_tcc_2 == 255][,1]
mask_del <- c(as.matrix(myraster$delaware_tcc_2 == 255))
coords_del <- crds(myraster)[mask_del,]

mydf_del <- data.frame(tcc = tcc_del, x = coords_del[,1],
                       y = coords_del[,2])
head(mydf_del)

# saveRDS(mydf_del, file = paste0(str_folder, "/tcc_delaware.rds"))
mydf_del <- readRDS(paste0(str_folder, "/tcc_delaware.rds"))

mygrid <- st_make_grid(st_bbox(shp), cellsize = 5000, square = F)
mygrid_fia <- mygrid[shp]
sample_fia <- st_sample(mygrid_fia, size = c(1,1))[shp]

plot(mygrid)
plot(st_geometry(shp), add = T)
plot(mygrid_fia, col = "#ff000088", add = T)
plot(sample_fia, cex = 0.3, add = T)

samples_coords <- st_coordinates(sample_fia)
samples_tcc <- exact_extract(myraster, st_buffer(sample_fia, 30),
                             "mean")$mean.delaware_tcc_1

mydf_del_samples <- data.frame(tcc = samples_tcc, x = samples_coords[,1],
                           y = samples_coords[,2])
ggplot(mydf_del_samples, aes(x = x, y = y, color = tcc)) + geom_point()

head(mydf_del_samples)

# saveRDS(mydf_del_samples, file = paste0(str_folder, "/tcc_del_samples.rds"))
mydf_del_samples <- readRDS(paste0(str_folder, "/tcc_del_samples.rds"))
 