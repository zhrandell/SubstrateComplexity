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

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## loads netCDF file
nc_data <- nc_open("T2015294191500.L2_LAC_OC.nc")


## saves text file with metadata 
{sink("T2015294191500.L2_LAC_OC_metadata.txt")
  print(nc_data)
  sink()}


## call variables 
## chlorophyll a data 
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
dim(chlor)


## pull lat and long data to anchor the chlorophyll observations 
slat <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/slat"))
slon <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/slon"))
elat <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/elat"))
elon <- as.data.frame(ncvar_get(nc_data, "scan_line_attributes/elon"))

dim(slat)
dim(slon)
dim(elat)
dim(elon)


## close nc_data file 
nc_close(nc_data)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=6,w=10, record=TRUE)

cuts = c(seq(0,2, length=20))
pal <- colorRampPalette(c("green","red"))


## this plots the entire "granule", a large region offshore of CA
p1 <- raster(t(chlor), 
            xmn=min(slon), xmx=max(elon), 
            ymn=min(slat), ymx=max(elat), 
            crs=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs+ towgs84=0,0,0"))

plot(p1, breaks=cuts, col=pal(20))


## specify rectangular region 
box1 <- as(extent(-119.9199, -119.0736, 32.9596, 33.5081), 'SpatialPolygons')
crs(box1) <- "+proj=longlat +datum=WGS84 +no_defs"


## apply cutout to pa 
p2 <- crop(p1, box1)
plot(p2, breaks=cuts, col=pal(20))


 
#download.file("http://www.naturalearthdata.com/http//www.naturalearthdata.com/download/10m/physical/ne_10m_coastline.zip", 
#              destfile = 'coastlines.zip')
#unzip(zipfile = "coastlines.zip", 
#      exdir = 'ne-coastlines-10m')

## add map of coastlines 
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
crs(coastlines) <- "+proj=longlat +datum=WGS84 +no_defs"
plot(coastlines)


## apply crop to coastline map to produce SNI image
SNI <- crop(coastlines, box1)
plot(SNI)


## plot both chlorophyll a and SNI 
plot(p2)
plot(SNI, add=TRUE)

## END script (for now) ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
