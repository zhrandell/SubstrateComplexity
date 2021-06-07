## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## test how to visualize / work with chlorophyll a data ~~~~~~~~~~~~~~~~~~~~~ ##
## Data obtained from NASA's MODIS Terra sensor, publicly available ~~~~~~~~~ ##
## Level 2 Chlorophyll data (not spatio-temporally averaged) ~~~~~~~~~~~~~~~~ ##
## https://oceancolor.gsfc.nasa.gov/ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## modified June 7th, 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library("raster")
library("ncdf4")
library("rgdal")
library("ggplot2")
library("RColorBrewer")
library("sp")  # classes for spatial data
library("rasterVis")  # raster visualisation
library("maptools")
library("rgeos")
library("dismo")
library("sf")
library("remotes")
library("gdalUtils")
library("rgeos")

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll")

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))

graphics.off()
windows(h=8,w=8, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## open and plot netCDF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
nc_data <- nc_open("grain3.nc")


## saves text file with metadata 
{sink("grain3_metadata.txt")
  print(nc_data)
  sink()}


## extract chlorophyll, latitude, and longitude information
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
lat <- ncvar_get(nc_data, "navigation_data/latitude")
long <- ncvar_get(nc_data, "navigation_data/longitude")


## close netCDF file 
nc_close(nc_data)


## open coastline data
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## key step I was previously missing: lat and long are current not in a raster 
## format, and should instead be treated as vectors. 
xy <- cbind(as.vector(long), as.vector(lat))


## visualize these vectors, where the lat long coordinates plot as the skewed 
## image that the previous raster format was unable to achieve
i <- seq(1, nrow(xy), 500)
plot(xy[i,])


## create RasterLayer out of lat / long coordinates 
## the stated resolution of r when resolution is undefined 
## (verify the above with r <- raster(extent(xy)), then run "r")
natural_res <- c(2.84775, 2.12091)


## scale the natural resolution by res
res <- 115
new_res <- natural_res / res


## creater raster layer
r <- raster(extent(xy), res=new_res)


## transform chlor values into Raster Layer 
x <- rasterize(xy, r, as.vector(chlor), 
               fun=mean)


## broader Channel Islands and southern CA bight
box1 <- as(extent(-121.0, -118.0, 32.4, 35.0), 'SpatialPolygons')


## zoom into SNI
box2 <- as(extent(-119.9199, -119.0736, 32.9596, 33.5081), 'SpatialPolygons')


## apply cutout to chlor, lat long raster & coastline data  
p1 <- crop(x, box2)
SNI <- crop(coastlines, box2)
plot(p1, zlim=c(0,10))
lines(SNI)

## MODIS implicit CRS -- explicitly set
MODIS_crs <- "+proj=longlat +datum=WGS84 +no_defs" 
crs(x) <- MODIS_crs


## check crs
#st_crs(x)
#st_crs(SNI)


## plot 
plot(p1, zlim=c(0,10))
lines(SNI)
abline(v=-119.4992, h=33.2465, col="blue")
## END additional data processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## create buffer areas of regions around SNI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create 2km buffer 
two_km_buffer <- raster::buffer(SNI, width=.03)
plot(two_km_buffer, col="cyan")
lines(SNI)


## use raster::erase to only retain outside of island
## need to convert SpatialLinesDataFrame --> sf --> SpatialPolygonDataFrame
SNI_sf <- st_as_sf(SNI)
SNI_poly <- st_polygonize(SNI_sf)
SNI_spatialpoly <- as(SNI_poly, "Spatial") 
buff <- raster::erase(two_km_buffer, SNI_spatialpoly)
plot(buff, col="#03A89E88")


## plot Cholorphyll data and 2km buffer
plot(p1, zlim=c(0,10), main="3km buffer to extract Chlorophyll values from around the outside of SNI")
plot(buff, col="#03A89E50", add=TRUE)
## END buffer creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Extract and visualize Chlorophyll data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## use spatial polygon specified by coordinates to extract data 
d1 <- raster::extract(p1, buff)
d1 <- na.omit(as.data.frame(d1))
names(d1)[1]<-"dat"


## plot
v1 <- ggplot(d1, aes(dat)) + my.theme +
  geom_histogram(binwidth=.25, color="black", fill="#388E8E") +
  xlab("Chlorophyll_a mg / m^-3") 
print(v1)
## END data extract ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





