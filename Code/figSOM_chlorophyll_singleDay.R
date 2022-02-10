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

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data/Chlorophyll")

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


## extract chlorophyll, latitude, and longitude information
chlor <- ncvar_get(nc_data, "geophysical_data/chlor_a")
lat <- ncvar_get(nc_data, "navigation_data/latitude")
long <- ncvar_get(nc_data, "navigation_data/longitude")
year <- ncvar_get(nc_data, "scan_line_attributes/year")
day <- ncvar_get(nc_data, "scan_line_attributes/day")


## close netCDF file 
nc_close(nc_data)


## open coastline data
setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data/Chlorophyll")
coastlines <- readOGR("ne-coastlines-10m/ne_10m_coastline.shp")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## lat and long are current not in a raster format, and should instead be treated as vectors. 
xy <- cbind(as.vector(long), as.vector(lat))


## visualize these vectors, where the lat long coordinates plot
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
x <- rasterize(xy, r, as.vector(chlor), fun=mean)


## narrow around SNI 
box2 <- as(extent(-119.8, -119.2, 33.05, 33.45), 'SpatialPolygons')


## apply cutout: zoom into SNI
p2 <- crop(x, box2)
SNI2 <- crop(coastlines, box2)


## plot SNI 
plot(p2, zlim=c(0,90))
lines(SNI2)
scalebar(15, xy=c(-119.4, 33.09), type='bar', divs=3, below="km")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





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


## extract the three new regions 
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





## Extract and visualize Chlorophyll data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## use spatial polygon specified by coordinates to extract data 
n_chlor <- as.data.frame(raster::extract(p2, n))
sw_chlor <- as.data.frame(raster::extract(p2, sw))
se_chlor <- as.data.frame(raster::extract(p2, se))


## site column 
n_chlor$region <- "North"
sw_chlor$region <- "SouthWest"
se_chlor$region <- "SouthEast"


## rename column with chlorophyll values
names(n_chlor)[1]<-"dat"
names(sw_chlor)[1]<-"dat" 
names(se_chlor)[1]<-"dat"


## single data frame with all chlorophyll data
chlorophyll_dat <- rbind(n_chlor, sw_chlor, se_chlor)


## reorder factors into desired order for plot 
chlorophyll_dat$region <- factor(chlorophyll_dat$region, levels=c("North", "SouthWest", "SouthEast"))
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



## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
