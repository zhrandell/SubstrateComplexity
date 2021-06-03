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

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
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
st_crs(coastlines)


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



## <- st_transform(coastlines, "+4326")
## <- projectRaster(coastline, crs = "+init=EPSG:4326")



## specify rectangular region right around SNI 
box1 <- as(extent(-119.9199, -119.0736, 32.9596, 33.5081), 'SpatialPolygons')
p2 <- crop(p1, box1)
extent(p2)
extent(p1)


plot(p2)
plot(coastlines, add=TRUE)
crs(coastlines)


## USE THIS METHOD FOR EXTRACTED DATA FROM RASTERS 
## create new extent (box / polygon to crop raster)
new.extent <- extent(-119.9199, -119.0736, 32.9596, 33.5081)
class(new.extent)


## crop raster 
p2 <- crop(x=p1, y=new.extent)
plot(p2)


## compare full and cropped raster 
plot(p1)
plot(p2, add=TRUE)


## use spatial polygon specified by coordinates to extract data 
d1 <- raster::extract(p1, new.extent)
d1 <- as.data.frame(d1)


## plot
v1 <- ggplot(d1, aes(d1)) + my.theme +
  geom_histogram(binwidth = .0025, color="black", fill="#388E8E") +
  xlab("Chlorophyll_a mg / m^-3") 
print(v1)





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
