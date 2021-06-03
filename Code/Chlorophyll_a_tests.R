## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## tests how to visualize / work with chlorophyll a data ~~~~~~~~~~~~~~~~~~~~ ##
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

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## loads netCDF file
#nc_data <- nc_open("T2015294191500.L2_LAC_OC.nc")
#nc_data <- nc_open("grain1.nc")
nc_data <- nc_open("grain2.nc")


## saves text file with metadata 
{sink("grain2_metadata.txt")
  print(nc_data)
  sink()}


## call variables 
## chlorophyll a data 
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")


lat <- as.data.frame(ncvar_get(nc_data, "navigation_data/latitude"))
long <- as.data.frame(ncvar_get(nc_data, "navigation_data/longitude")) 


## pull lat and long data to anchor the chlorophyll observations 
#slat <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/slat"))
#slon <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/slon"))
#elat <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/elat"))
#elon <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/elon"))

#dim(slat)
#dim(slon)
#dim(elat)
#dim(elon)


## close nc_data file 
nc_close(nc_data)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=7,w=8, record=TRUE)


## add map of coastlines 
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
#projection(coastlines) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(coastlines)
coastlines



## this plots the entire "granule", a large region offshore of CA
p1 <- raster(t(chlor), 
            xmn=min(long), xmx=max(long), 
            ymn=min(lat), ymx=max(lat))


## coordinate system information 
#crs=crs("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"))
## set crs
#projection(p1) <- "+proj=longlat +datum=WGS84 +no_defs"
crs(p1)


## plot entire granule 
plot(p1, zlim=c(0,1), main="entire granule (southern CA in far upper right corner)")
plot(coastlines, add=TRUE)



## specify rectangular region right around SNI 
box1 <- as(extent(-119.9199, -119.0736, 32.9596, 33.5081), 'SpatialPolygons')
p2 <- crop(p1, box1)

plot(p2)
plot(coastlines, add=TRUE)
crs(coastlines)

## specify broader region offshore of southern CA
box2 <- as(extent(-121.0, -117.0, 32.4, 38.0), 'SpatialPolygons')
#projection(box2) <- "+proj=longlat +datum=WGS84 +no_defs"

## exact box from NASA site
box3 <- as(extent(-120.6162, -118.2973, 32.7248, 34.3512), 'SpatialPolygons')




## apply cutout to pa 
p3 <- crop(p1, box3)
plot(p3, zlim=c(0,1), main="chlor_a data mismatch in coordinate system")
plot(coastlines, add=TRUE)



p2
crs(coastlines)
crs(box2)


## END script (for now) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
