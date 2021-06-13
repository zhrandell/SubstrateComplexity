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

#save.image(file='MODIS_Terra_7685files.RData')
#quit(save='no')

## load previous spatial analysis 
#load('MODIS_Aqua_1626files.RData')
#load('MODIS_Aqua_7685files.RData')


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
windows(h=5,w=8, record=TRUE)
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

## 1 degree latitude --> km  
km <- 1 / 110.574

## create buffer areas of regions around SNI ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
buffer <- raster::buffer(SNI, width=3*km)

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
datasets <- list.files('D:/OneDrive/Active_Projects/Substrate_Complexity/Data/Chlorophyll/MODIS_Aqua/testRun', 
                       pattern="*.nc", full.names=TRUE, include.dirs=FALSE)
## an alternative
#filelist <- paste("automate/", dir("automate"), sep="")


## set up empty numeric vectors to hold Chlorophyll data 
NW <- numeric()
NE <- numeric()
SE <- numeric()
SW <- numeric()
year <- numeric()
day <- numeric()

## number of NETcdf files to loop through
N <- length(datasets)

## for loop 
for (i in 1:N){
  
    ## open NETcdf file, extract relevant information, and close the file
    nc_data <- nc_open(datasets[i])
    chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
    lat <- ncvar_get(nc_data, "navigation_data/latitude")
    long <- ncvar_get(nc_data, "navigation_data/longitude")
    yr <- ncvar_get(nc_data, "scan_line_attributes/year")
    dy <- ncvar_get(nc_data, "scan_line_attributes/day")
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
    
    year[i] <- list(rep(yr[1], length(NW[[i]])))
    day[i] <- list(rep(dy[1], length(NW[[i]])))
    
}
## END for loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


## post-loop data processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## expand list of lists into single list and convert to data frame
NW.list <- na.omit(as.data.frame(unlist(NW)))
NE.list <- na.omit(as.data.frame(unlist(NE)))
SE.list <- na.omit(as.data.frame(unlist(SE)))
SW.list <- na.omit(as.data.frame(unlist(SW)))

## Day and Year extraction and data processing grayed out for now 
#NW.list <- as.data.frame(unlist(NW))
#NW.list$logChlor <- log(NW.list[,1])
#NW.day <- as.data.frame(unlist(day))
#NW.year <- as.data.frame(unlist(year))
#new_NW <- na.omit(cbind(NW.list, NW.day, NW.year))
#names(new_NW)[1]<-"chlor"
#names(new_NW)[3]<-"day" 
#names(new_NW)[4]<-"year" 
#new_NW$region <- "NorthWest"
#head(new_NW)

## sample sizes
len.NW <- length(NW.list[,1])
len.NE <- length(NE.list[,1])
len.SW <- length(SW.list[,1])
len.SE <- length(SE.list[,1])

## minimum sample size
sampleSize <- min(len.NW, len.NE, len.SW, len.SE)

## randomly sample from data.frames to ensure equal sample size
NorthWest <- as.data.frame(sample(NW.list[,1], sampleSize, replace = FALSE))
NorthEast <- as.data.frame(sample(NE.list[,1], sampleSize, replace = FALSE))
SouthWest <- as.data.frame(sample(SW.list[,1], sampleSize, replace = FALSE))
SouthEast <- as.data.frame(sample(SE.list[,1], sampleSize, replace = FALSE))

## rename column with chlorophyll values
names(NorthWest)[1]<-"chlor" 
names(NorthEast)[1]<-"chlor"
names(SouthWest)[1]<-"chlor"
names(SouthEast)[1]<-"chlor"

## log scale
NorthWest$logChlor <- log(NorthWest$chlor)
NorthEast$logChlor <- log(NorthEast$chlor)
SouthWest$logChlor <- log(SouthWest$chlor)
SouthEast$logChlor <- log(SouthEast$chlor)

## add region column 
NorthWest$region <- "NorthWest"
NorthEast$region <- "NorthEast"
SouthWest$region <- "SouthWest"
SouthEast$region <- "SouthEast"


## bind into single data frame
dat <- rbind(NorthWest, NorthEast, SouthWest, SouthEast)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot on normal and log scale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## set x and y lim here 
## single kernal density plot
my.pal <- c("#476A34", "#FF6103", "#B22222", "#008080")

## arrange all plots in a single pane
graphics.off()
windows(h=5,w=8, record=TRUE)


## reorder factors into desired order for plot 
dat$region <- factor(dat$region, levels=c("NorthWest", "NorthEast","SouthWest", "SouthEast"))



log.KD <- ggplot(dat, aes(logChlor, fill=region)) + my.theme +
  geom_density(alpha=.4) +  xlab("log Chlorophyll mg / m^-3")
print(log.KD)



KD <- ggplot(dat, aes(chlor, fill=region)) + my.theme +
  geom_density(alpha=.4) +  xlab("log Chlorophyll mg / m^-3") 

print(KD)
## END kernal density and velocity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


range(dat$chlor)


## Inverse CDF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## calcuate empirical cumulative density function and extract and sort values
NW_ecdf <- data.frame(x=unique(NorthWest$logChlor), 
                      y=ecdf(NorthWest$logChlor)(unique(NorthWest$logChlor))*length(NorthWest$logChlor))
NE_ecdf <- data.frame(x=unique(NorthEast$logChlor), 
                      y=ecdf(NorthEast$logChlor)(unique(NorthEast$logChlor))*length(NorthEast$logChlor))

SW_ecdf <- data.frame(x=unique(SouthWest$logChlor), 
                      y=ecdf(SouthWest$logChlor)(unique(SouthWest$logChlor))*length(SouthWest$logChlor))
SE_ecdf <- data.frame(x=unique(SouthEast$logChlor), 
                      y=ecdf(SouthEast$logChlor)(unique(SouthEast$logChlor))*length(SouthEast$logChlor))



## rescale extracted cdf values to 0-1 scale
NW_ecdf$y <- scale(NW_ecdf$y, center=min(NW_ecdf$y), scale=diff(range(NW_ecdf$y)))
NE_ecdf$y <- scale(NE_ecdf$y, center=min(NE_ecdf$y), scale=diff(range(NE_ecdf$y)))
SW_ecdf$y <- scale(SW_ecdf$y, center=min(SW_ecdf$y), scale=diff(range(SW_ecdf$y)))
SE_ecdf$y <- scale(SE_ecdf$y, center=min(SE_ecdf$y), scale=diff(range(SE_ecdf$y)))


## take the inverse of a cdf, such that the p(x) > or = log wave event
NW_ecdf$inv_y <- ((NW_ecdf$y - max(NW_ecdf$y)) * (-1))
NE_ecdf$inv_y <- ((NE_ecdf$y - max(NE_ecdf$y)) * (-1))
SW_ecdf$inv_y <- ((SW_ecdf$y - max(SW_ecdf$y)) * (-1))
SE_ecdf$inv_y <- ((SE_ecdf$y - max(SE_ecdf$y)) * (-1))


## add site name to new data frames
NW_ecdf$region <- "NorthWest"
NE_ecdf$region <- "NorthEast"
SW_ecdf$region <- "SouthWest"
SE_ecdf$region <- "SouthEast"


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat_ecdf <- rbind(NW_ecdf, NE_ecdf, SW_ecdf, SE_ecdf)
dat_ecdf$region <- factor(dat_ecdf$region, levels=c("NorthWest", "NorthEast", "SouthWest", "SouthEast"))


## plot all 
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=region)) + my.theme +
  geom_line(lwd=1, alpha=.7) +
  #scale_color_manual(values=pal_sites) +
  xlab("log Chlorophyll mg / m^-3") + ylab("inverse empirical CDF")  

print(p1)
## End script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







