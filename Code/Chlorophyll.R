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
library("lwgeom")

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll/automate")

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))

graphics.off()
windows(h=6.5,w=8, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## open and plot netCDF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
nc_data <- nc_open("grain4.nc")

## saves text file with metadata 
#{sink("grain2_metadata.txt")
#  print(nc_data)
#  sink()}

## extract chlorophyll, latitude, and longitude information
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
lat <- ncvar_get(nc_data, "navigation_data/latitude")
long <- ncvar_get(nc_data, "navigation_data/longitude")

## close netCDF file 
nc_close(nc_data)

## open coastline data
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll")
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
res <- 120
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

## narrow around SNI 
box3 <- as(extent(-119.8, -119.2, 33.05, 33.45), 'SpatialPolygons')

## apply cutout to chlor, lat long raster & coastline data  
p1 <- crop(x, box3)
SNI <- crop(coastlines, box3)

plot(x, zlim=c(0,10))
lines(SNI)

plot(p1)

## MODIS implicit CRS -- explicitly set
#MODIS_crs <- "+proj=longlat +datum=WGS84 +no_defs" 
#crs(x) <- MODIS_crs


## plot 
plot(p1)
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
plot(p1,  
     main="3km buffer to extract Chlorophyll values from around the outside of SNI")
plot(buff, col="#03A89E50", add=TRUE)
## END buffer creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## CREATE LINES TO PARTITION BUFFER AND CREATE FOUR REGIONS ~~~~~~~~~~~~~~~~~~~~
## coordinates centered around SNI lat / long 
#l1 <-  rbind(c(-119.4992, 33.1), c(-119.4992, 33.4))
#l2 <- rbind(c(-119.8, 33.2468), c(-119.2, 33.2468))


## coordinates manually determined lines based on exposure
l1 <-  rbind(c(-119.312897, 33.18), c(-119.695722, 33.32))
l2 <- rbind(c(-119.441736, 33.38125), c(-119.570142, 33.125417))

## coordintes --> Formal class Line object 
S1 <- Line(l1)
S2 <- Line(l2)
S3 <- Lines(list(S1, S2), ID="a")

## Formal class Line --> SpatialLines 
quadLines <- SpatialLines(list(S3))

## check lines 
plot(p1)
lines(SNI)
plot(buff, add=T)
lines(quadLines)

## remove coordinate system from buffer layer to facilitate next steps
crs(buff) <- NA
crs(quadLines) <- NA

## intersect buffer with quadLines, buffer the intersection, and split the regions
intersect1 <- gIntersection(buff, quadLines)
intersect2 <- gBuffer(intersect1, width=0.0001) 
regions <- gDifference(buff, intersect2)

## extract the four new regions 
sw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[1]]), "1")))
nw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[2]]), "2")))
ne <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[3]]), "3")))
se <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[4]]), "4")))

## plot and check
plot(p1, main="3km buffer split into four regions")
lines(SNI)
plot(sw, add=T, col="#B2222250")
plot(nw, add=T, col="#476A3450")
plot(ne, add=T, col="#FF610350") 
plot(se, add=T, col="#00808050")
## END buffer regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Extract and visualize Chlorophyll data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## use spatial polygon specified by coordinates to extract data 
sw_chlor <- as.data.frame(raster::extract(p1, sw))
nw_chlor <- as.data.frame(raster::extract(p1, nw))
ne_chlor <- as.data.frame(raster::extract(p1, ne))
se_chlor <- as.data.frame(raster::extract(p1, se))

## site column 
sw_chlor$region <- "SouthWest"
nw_chlor$region <- "NorthWest"
ne_chlor$region <- "NorthEast"
se_chlor$region <- "SouthEast"

## rename column with chlorophyll values
names(sw_chlor)[1]<-"dat" 
names(nw_chlor)[1]<-"dat"
names(ne_chlor)[1]<-"dat"
names(se_chlor)[1]<-"dat"

## single data frame with all chlorophyll data
chlorophyll_dat <- rbind(sw_chlor, nw_chlor, ne_chlor, se_chlor)

## reorder factors into desired order for plot 
chlorophyll_dat$region <- factor(chlorophyll_dat$region, levels=c("NorthWest", "NorthEast",
                                                                  "SouthWest", "SouthEast"))
## END chlorophyll data extract ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~           

                                                       


                                                                                                           
## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
graphics.off()
windows(h=6,w=10, record=TRUE)

## custom color pallet
my.pal <- c("#476A34","#FF6103",
            "#B22222", "#008080")

## frequency histograms 
v1 <- ggplot(chlorophyll_dat, aes(dat, fill=region)) + geom_histogram(binwidth=.5) +
  xlab("Chlorophyll_a mg / m^-3") + my.theme + facet_wrap(~ region) +
  scale_fill_manual(values = my.pal)
print(v1)
## END data extract ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##

