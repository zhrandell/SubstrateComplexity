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

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data/Chlorophyll")


## load previous spatial analysis 
#load('MODIS_Aqua_1626files.RData')
load('MODIS_Terra_7685files.RData')


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





## crop data and create buffer part 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Select one -- 1: broader Channel Island, 2: SNI medium view, 3: SNI zoomed in ht
#box1 <- as(extent(-121.0, -118.0, 32.4, 35.0), 'SpatialPolygons')
box2 <- as(extent(-119.8, -119.2, 33.05, 33.45), 'SpatialPolygons')

## open coastline data
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")

## crop coastlines to create buffer
SNI <- crop(coastlines, box2)

natural_res <- c(2.84775, 2.12091)
res <- 120
new_res <- natural_res / res
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## use pts to create 3km buffer around sites ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac, WestEnd, EastDutch, Daytona
x_lon <- c(-119.48681, -119.57419, -119.48407, -119.44412)
y_lat <- c(33.27354, 33.24772, 33.21598, 33.21687)
latlong <- cbind(x_lon, y_lat)
pts <- SpatialPoints(latlong)
proj4string(pts) <- CRS("+init=epsg:4326")

## create 3km buffer around pts and check via plot 
pt_buffer <- buffer(pts, 3000)
plot(pt_buffer, col="cyan")
lines(SNI)
points(pts)

## intersect and cutout island points
## remove coordinate system from buffer layer to facilitate next steps
crs(pt_buffer) <- NA
crs(SNI) <- NA
intersect <- gIntersection(pt_buffer, SNI)
intersect2 <- gBuffer(intersect, width=0.00001)
regions <- gDifference(pt_buffer, intersect2)

## extract the four new regions 
sw <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[3]]), "3")))
n <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[4]]), "4")))
se <- SpatialPolygons(list(Polygons(list(regions@polygons[[1]]@Polygons[[7]]), "7")))


## plot and check
plot(pt_buffer, border="white")
plot(sw, add=T, col="#C77CFF50")
plot(n, add=T, col="#7CAE0050")
plot(se, add=T, col="#F8766D50")
lines(SNI)
points(pts, pch=21, col="black", bg="red", cex=1)
## end pt buffer creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Loop through all NETcdf files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## create file list
datasets <- list.files('D:/OneDrive/Active_Projects/SubstrateComplexity/Data/Chlorophyll/MODIS_Terra/files', 
                       pattern="*.nc", 
                       full.names=TRUE, 
                       include.dirs=FALSE)


## set up empty numeric vectors to hold Chlorophyll data 
N <- numeric()
SW <- numeric()
SE <- numeric()
year <- numeric()
day <- numeric()

## number of NETcdf files to loop through
len <- length(datasets)

## loop through all NETcdf files 
for (i in 1:len){
  
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
    p1 <- crop(x, box2)
    
    ## chlorophyll data extraction
    N[i] <- raster::extract(p1, n)
    SW[i] <- raster::extract(p1, sw)
    SE[i] <- raster::extract(p1, se)
    
    year[i] <- list(rep(yr[1], length(N[[i]])))
    day[i] <- list(rep(dy[1], length(N[[i]])))
    
}
## END NETcdf for loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## post-loop data processing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## expand list of lists into single list and convert to data frame
N.list <- na.omit(as.data.frame(unlist(N)))
SE.list <- na.omit(as.data.frame(unlist(SE)))
SW.list <- na.omit(as.data.frame(unlist(SW)))


## sample sizes
len.N <- length(N.list[,1])
len.SW <- length(SW.list[,1])
len.SE <- length(SE.list[,1])


## minimum sample size
sampleSize <- min(len.N, len.SW, len.SE)


## randomly sample from data.frames to ensure equal sample size
North <- as.data.frame(sample(N.list[,1], sampleSize, replace = FALSE))
SouthWest <- as.data.frame(sample(SW.list[,1], sampleSize, replace = FALSE))
SouthEast <- as.data.frame(sample(SE.list[,1], sampleSize, replace = FALSE))


## rename column with chlorophyll values
names(North)[1]<-"chlor" 
names(SouthWest)[1]<-"chlor"
names(SouthEast)[1]<-"chlor"


## log scale
North$logChlor <- log10(North$chlor)
SouthWest$logChlor <- log10(SouthWest$chlor)
SouthEast$logChlor <- log10(SouthEast$chlor)


## add region column 
North$region <- "North"
SouthWest$region <- "SouthWest"
SouthEast$region <- "SouthEast"


## bind into single data frame
dat <- rbind(North, SouthWest, SouthEast)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot on normal and log scale ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## single kernal density plot
my.pal <- c("#7CAE0050", "#C77CFF50", "#F8766D50")

## arrange all plots in a single pane
graphics.off()
windows(h=5,w=8, record=TRUE)

## reorder factors into desired order for plot 
dat$region <- factor(dat$region, levels=c("North", "SouthWest", "SouthEast"))

## plot on the log scale
log.KD <- ggplot(dat, aes(logChlor, fill=region)) + my.theme +
  scale_fill_manual(values=my.pal) +
  geom_density(alpha=.4) +  xlab("Chlorophyll mg / m^-3") +
  scale_x_continuous(labels=c("0.001", "0.01", "0.1", "1", "10", "100", "1000"))
print(log.KD)

## plot on the natural scale 
KD <- ggplot(dat, aes(chlor, fill=region)) + my.theme +
  geom_density(alpha=.4) +  xlab("log Chlorophyll mg / m^-3") 

print(KD)
## END kernal density and velocity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
N_num <- North$logChlor
SW_num <- SouthWest$logChlor
SE_num <- SouthEast$logChlor

## perform KS tests 
ks.test(N_num, SW_num)
ks.test(N_num, SE_num)
ks.test(SW_num, SE_num)
## End KS tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Inverse CDF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## t1 = calculate empirical cumulative density function values
## t2 = rescale extracted cdf values to 0-1 scale
## t3 = take the inverse of a cdf, such that the p(x) > or = log wave event
ecdf.f <- function(dat, site.name){
  t1 <- data.frame(x=unique(dat$logChlor), y=ecdf(dat$logChlor)(unique(dat$logChlor))*length(dat$logChlor))
  t2 <- scale(t1$y, center=min(t1$y), scale=diff(range(t1$y)))
  t3 <- ((t2 - max(t2)) * (-1))
  t4 <- cbind(t1, t2, t3)
  names(t4)[1:4]<-c("x", "unscaled_y", "y", "inv_y")
  t4$Site <- site.name 
  out <- return(t4)
}


## perform calculations
N_ecdf <- ecdf.f(NF, "North")
SW_ecdf <- ecdf.f(WE, "SouthWest")
SE_ecdf <- ecdf.f(Day, "SouthEast")


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat_ecdf <- rbind(N_ecdf, SW_ecdf, SE_ecdf)
dat_ecdf$Site <- factor(dat_ecdf$region, levels=c("North", "SouthWest", "SouthEast"))
## END ecdf calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## inverse ecdf plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=region)) + my.theme +
  geom_line(lwd=1, alpha=.7) + 
  scale_color_manual(values=my.pal) +
  xlab("Chlorophyll mg / m^-3") + ylab("inverse empirical CDF") +
  scale_x_continuous(n.breaks=7, labels=c("0.001", "0.01", "0.1", "1", "10", "100", "1000"))

print(p1)
## End of plots~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
