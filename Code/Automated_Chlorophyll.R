## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## test how to visualize / work with chlorophyll a data ~~~~~~~~~~~~~~~~~~~~~ ##
## Data obtained from NASA's MODIS Terra sensor, publicly available ~~~~~~~~~ ##
## Level 2 Chlorophyll data (not spatio-temporally averaged) ~~~~~~~~~~~~~~~~ ##
## https://oceancolor.gsfc.nasa.gov/ ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## modified June 7th, 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## start up ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library("ggpubr")
library("raster")
library("ncdf4")
library("ncdf4.helpers")
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
library("RNetCDF")
library("rlist")
library("reshape")
library("egg")

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
windows(h=6.5,w=8, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




#####  

## FIXED DEFINITIONS AND OBJECTS USED FOR ALL FILES ~~~~~~~~~~~~~~~~~~~~~~~~

## CROP DATA AND CREATE BUFFER PART 1~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Select one -- 1: broader Channel Island, 2: SNI medium view, 3: SNI zoomed in ht
#box1 <- as(extent(-121.0, -118.0, 32.4, 35.0), 'SpatialPolygons')
#box2 <- as(extent(-119.9199, -119.0736, 32.9596, 33.5081), 'SpatialPolygons')
box3 <- as(extent(-119.8, -119.2, 33.05, 33.45), 'SpatialPolygons')

## open coastline data
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")

## crop coastlines to create buffer
SNI <- crop(coastlines, box3)

## create buffer areas of regions around SNI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
buffer <- raster::buffer(SNI, width=.03)

## use raster::erase to only retain outside of island
SNI_sf <- st_as_sf(SNI)
SNI_poly <- st_polygonize(SNI_sf)
SNI_spatialpoly <- as(SNI_poly, "Spatial") 
buff <- raster::erase(buffer, SNI_spatialpoly)
plot(buff, col="cyan")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## CREATE LINES TO DELINIATE REGIONS WITHIN BUFFER ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## coordinates manually determined 
l1 <-  rbind(c(-119.312897, 33.18), c(-119.695722, 33.32))
l2 <- rbind(c(-119.441736, 33.38125), c(-119.570142, 33.125417))

## coordinates --> Formal class Line object 
S1 <- Line(l1)
S2 <- Line(l2)
S3 <- Lines(list(S1, S2), ID="a")

## Formal class Line --> SpatialLines 
quadLines <- SpatialLines(list(S3))

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
plot(buff)
plot(sw, add=T, col="#B2222250")
plot(nw, add=T, col="#476A3450")
plot(ne, add=T, col="#FF610350") 
plot(se, add=T, col="#00808050")
lines(quadLines)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## FINAL DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## set implicit crs for Chlor raster layer
MODIS_crs <- "+proj=longlat +datum=WGS84 +no_defs" 

## create RasterLayer out of lat / long coordinates 
natural_res <- c(2.84775, 2.12091)

## scale the natural resolution by res
res <- 120
new_res <- natural_res / res
## END UNIVERSAL OBJECTS AND DEFINITIONS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#####






## Loop through all NETcdf files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

## create file list
datasets <- list.files('D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll/automate', 
                       pattern="*.nc", full.names=TRUE, include.dirs=FALSE)
## an alternative
#filelist <- paste("automate/", dir("automate"), sep="")




#out<-dim(0)




## set up empty numeric vectors to hold Chlorophyll data 
NW <- numeric()
NE <- numeric()
SE <- numeric()
SW <- numeric()

## number of NETcdf files to loop through
N <- length(datasets)
for (i in 1:N){
  
    ## open NETcdf file, extract relevant information, and close the file
    nc_data <- nc_open(datasets[i])
    chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
    lat <- ncvar_get(nc_data, "navigation_data/latitude")
    long <- ncvar_get(nc_data, "navigation_data/longitude")
    nc_close(nc_data)
    
    ## data processing
    xy <- cbind(as.vector(long), as.vector(lat))
    r <- raster(extent(xy), res=new_res)
    x <- rasterize(xy, r, as.vector(chlor), fun=mean)
    p1 <- crop(x, box3)
    
    ## chlorophyll data extraction
    NW[i] <- raster::extract(p1, nw)
    NE[i] <- raster::extract(p1, ne)
    SE[i] <- raster::extract(p1, se)
    SW[i] <- raster::extract(p1, sw)
}
    

## expand list of lists into single list and convert to data frame
NW.list <- as.data.frame(unlist(NW))
NE.list <- as.data.frame(unlist(NE))
SE.list <- as.data.frame(unlist(SE))
SW.list <- as.data.frame(unlist(SW))


## rename column with chlorophyll values
names(NW.list)[1]<-"dat" 
names(NE.list)[1]<-"dat"
names(SE.list)[1]<-"dat"
names(SW.list)[1]<-"dat"

## individual region plots 
plot.NW <- ggplot(NW.list, aes(dat)) + geom_histogram(binwidth=.15, fill="#476A34") + my.theme + xlim(0,15) + ylim(0,15) + ggtitle("NorthWest") + xlab("Chlorophyll_a mg / m^-3") + 
  theme(axis.title.x = element_blank(), axis.text.x = element_blank())
plot.NE <- ggplot(NE.list, aes(dat)) + geom_histogram(binwidth=.15, fill="#FF6103") + my.theme + xlim(0,15) + ylim(0,15) + ggtitle("NorthEast") + xlab("Chlorophyll_a mg / m^-3") +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank())
plot.SW <- ggplot(SW.list, aes(dat)) + geom_histogram(binwidth=.15, fill="#B22222") + my.theme + xlim(0,15) + ylim(0,15) + ggtitle("SouthWest") + xlab("Chlorophyll_a mg / m^-3") 
plot.SE <- ggplot(SE.list, aes(dat)) + geom_histogram(binwidth=.15, fill="#008080") + my.theme + xlim(0,15) + ylim(0,15) + ggtitle("SouthEast") + xlab("Chlorophyll_a mg / m^-3") +
  theme(axis.title.y = element_blank(), axis.text.y = element_blank())

## arrange all plots in a single pane
graphics.off()
windows(h=6.5,w=12, record=TRUE)

g1 <- ggarrange(tag_facet(plot.NW + facet_wrap(~"dat"), tag_pool="a"), 
                tag_facet(plot.NE + facet_wrap(~"dat"), tag_pool="b"),
                tag_facet(plot.SW + facet_wrap(~"dat"), tag_pool="c"),
                tag_facet(plot.SE + facet_wrap(~"dat"), tag_pool="d"),
                nrow=2)

print(g1)
## END kernal density and velocity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








