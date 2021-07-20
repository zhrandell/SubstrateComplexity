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
library("dismo")
library("sf")
library("remotes")
library("gdalUtils")
library("rgeos")
library("lwgeom")

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll")

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16), 
                 legend.text = element_text(size=14), 
                 legend.title = element_blank(), 
                 legend.position = c(0.85, 0.85))

graphics.off()
windows(h=5,w=6, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## open and plot netCDF file ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
nc_data <- nc_open("grain3.nc")

## saves text file with metadata 
#{sink("grain2_metadata.txt")
#  print(nc_data)
#  sink()}

## extract chlorophyll, latitude, and longitude information
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
lat <- ncvar_get(nc_data, "navigation_data/latitude")
long <- ncvar_get(nc_data, "navigation_data/longitude")
year <- ncvar_get(nc_data, "scan_line_attributes/year")
day <- ncvar_get(nc_data, "scan_line_attributes/day")




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
#box1 <- as(extent(-121.0, -118.0, 32.4, 35.0), 'SpatialPolygons')

## narrow around SNI 
box2 <- as(extent(-119.8, -119.2, 33.05, 33.45), 'SpatialPolygons')

## apply cutout: broader Channel Islands region  
#p1 <- crop(x, box1)
#SNI1 <- crop(coastlines, box1)

## apply cutout: zoom into SNI
p2 <- crop(x, box2)
SNI2 <- crop(coastlines, box2)

## plot SNI 
plot(p2, zlim=c(0,90))
lines(SNI2)
scalebar(15, xy=c(-119.4, 33.09), type='bar', divs=3, below="km")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## create buffer areas of regions around SNI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## 1 degree latitude --> km  
#km <- 1 / 110.574

## create 2km buffer 
#buffer <- raster::buffer(SNI2, width=3*km)
#plot(buffer, col="cyan")
#lines(SNI2)


#buffer2 <- buffer()

## use raster::erase to only retain outside of island
## need to convert SpatialLinesDataFrame --> sf --> SpatialPolygonDataFrame
#SNI_sf <- st_as_sf(SNI2)
#SNI_poly <- st_polygonize(SNI_sf)
#SNI_spatialpoly <- as(SNI_poly, "Spatial") 
#buff <- raster::erase(buffer, SNI_spatialpoly)
#plot(buff, col="#03A89E88")

## plot Cholorphyll data and 2km buffer
#plot(p2)
#plot(buff, col="#03A89E50", add=TRUE)
## END buffer creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## use site lat long to create raster-extraction regions ~~~~~~~~~~~~~~~~~~~~~~~

## NavFac, WestEnd, EastDutch, Daytona
x_lon <- c(-119.48681, -119.57419, -119.48407, -119.44412)
y_lat <- c(33.27354, 33.24772, 33.21598, 33.21687)
latlong <- cbind(x_lon, y_lat)
pts <- SpatialPoints(latlong)
proj4string(pts) <- CRS("+init=epsg:4326")

## create 3km buffer around pts and check via plot 
pt_buffer <- buffer(pts, 3000)
plot(p2)
plot(pt_buffer, col="cyan", add=T)
points(pts)
lines(SNI2)

## intersect and cutout island points
## remove coordinate system from buffer layer to facilitate next steps
crs(pt_buffer) <- NA
crs(SNI2) <- NA
intersect <- gIntersection(pt_buffer, SNI2)
intersect2 <- gBuffer(intersect, width=0.00001)
regions <- gDifference(pt_buffer, intersect2)

## extract the four new regions 
sw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[3]]), "3")))
n <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[4]]), "4")))
se <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[7]]), "7")))


## plot and check
pal <- colorRampPalette(c("white", "#0099CC", "yellow"))

plot(p2, col=pal(1000))
lines(SNI2)
plot(sw, add=T, col="#C77CFF50")
plot(n, add=T, col="#7CAE0050")
plot(se, add=T, col="#F8766D50")
points(pts, pch=21, col="black", bg="red", cex=1)
#scalebar(12, xy=c(-119.4, 33.09), type='bar', divs=5, lonlat = FALSE, label=c('0','6','12'), below="km")
compassRose(-119.35,33.385, cex=.75)
text(-119.48681, 33.27354, "NavFac", pos=3, offset=2)
text(-119.57419, 33.24772, "WestEnd", pos=2, offset=1.45)
text(-119.48547, 33.21652, "DutchHarbor", pos=1, offset=1.6)
text(-119.44412, 33.21687, "Daytona", pos=4, offset=1.65)

## end point creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## CREATE LINES TO PARTITION BUFFER AND CREATE FOUR REGIONS ~~~~~~~~~~~~~~~~~~~~
## coordinates centered around SNI lat / long 
#l1 <-  rbind(c(-119.4992, 33.1), c(-119.4992, 33.4))
#l2 <- rbind(c(-119.8, 33.2468), c(-119.2, 33.2468))


## coordinates manually determined lines based on exposure
#l1 <-  rbind(c(-119.312897, 33.18), c(-119.695722, 33.32))
#l2 <- rbind(c(-119.441736, 33.38125), c(-119.570142, 33.125417))

## coordintes --> Formal class Line object 
#S1 <- Line(l1)
#S2 <- Line(l2)
#S3 <- Lines(list(S1, S2), ID="a")

## Formal class Line --> SpatialLines 
#quadLines <- SpatialLines(list(S3))

## check lines 
#plot(p2)
#lines(SNI2)
#plot(buff, add=T)
#lines(quadLines)

## remove coordinate system from buffer layer to facilitate next steps
#crs(buff) <- NA
#crs(quadLines) <- NA

## intersect buffer with quadLines, buffer the intersection, and split the regions
#intersect1 <- gIntersection(buff, quadLines)
#intersect2 <- gBuffer(intersect1, width=0.0001) 
#regions <- gDifference(buff, intersect2)

## extract the four new regions 
#sw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[1]]), "1")))
#nw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[2]]), "2")))
#ne <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[3]]), "3")))
#se <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[4]]), "4")))

## plot and check
#plot(p2, main="3km buffer split into four regions")
#lines(SNI2)
#plot(sw, add=T, col="#C77CFF50")
#plot(nw, add=T, col="#7CAE0050")
#plot(ne, add=T, col="#F8766D50") 
#plot(se, add=T, col="#00BFC450")
#plot(buff, add=T)
## END buffer regions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## full map ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## GPS coordinates for SNI sites
#coords <- as.data.frame(t(array(c(
#  "NavFac", -119.48681, 33.27354,
#  "WestEnd", -119.57419, 33.24772,
#  "WestDutch", -119.48547, 33.21652,
#  "EastDutch", -119.48407, 33.21598,
#  "Daytona", -119.44412, 33.21687 
#), dim=c(3,5))))

## plot full map 
#graphics.off()
#windows(h=5,w=6, record=TRUE)
#legend.args = list(text="Chlorophyll mg / m^-3")
#plot(p2, zlim=c(0,60))
#lines(SNI2)
#plot(buffer, add=TRUE)
#plot(sw, add=T, col="#00BFC470")
#plot(nw, add=T, col="#F8766D70")
#plot(ne, add=T, col="#7CAE0070")
#plot(se, add=T, col="#C77CFF70")
#scalebar((12*km), xy=c(-119.4, 33.09), type='bar', divs=5, lonlat = FALSE, label=c('0','6','12'), below="km")
#points(x=coords[,2], y=coords[,3], pch=21, col="black", bg="red", cex=1)
#text(-119.48681, 33.27354, "NavFac", pos=3, offset=2)
#text(-119.57419, 33.24772, "WestEnd", pos=2, offset=1.45)
#text(-119.48547, 33.21652, "DutchHarbor", pos=1, offset=1.6)
#text(-119.44412, 33.21687, "Daytona", pos=4, offset=1.65)
#compassRose(-119.35,33.385, cex=.75)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Extract and visualize Chlorophyll data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## use spatial polygon specified by coordinates to extract data 
n_chlor <- as.data.frame(raster::extract(p1, n))
sw_chlor <- as.data.frame(raster::extract(p1, sw))
se_chlor <- as.data.frame(raster::extract(p1, se))
#se_chlor <- as.data.frame(raster::extract(p1, se))

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






