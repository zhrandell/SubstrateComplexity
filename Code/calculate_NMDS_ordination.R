#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## NON-METRIC MULTIDIMENSIONAL SCALING (NMDS) analyses for SNI subtidal data 
## Zachary Randell; code updated 10 November 2020
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## clear Global Environment 
rm(list = ls())

## load packages (if packages not installed: install.packages("librarynamehere"))
library(ggplot2)
library(dplyr)
library(vegan)
library(MASS)
library(fitdistrplus)

## close out window tab & open custom window
graphics.off()
windows(h=8,w=8, record=TRUE)

## set working directory
setwd("D:/Active_Projects/Substrate_Complexity/Code")

## load data frame
dat <- read.csv("SNI_subtidal_swath.csv", header = TRUE)

## load ordination (results of a previously run ordination, so one does not have to take large amounts of time to rerun analysis) 
load("ord_22March2020_noLOG.rda")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Data prep 
## set site as factor and rename 
dat$SITE <- factor(dat$SITE, levels=c("1", "2", "3", "4",
                                      "5", "6", "7"), 
                   labels=c("NavFac", "WestEnd Urchin", "WestEnd Kelp",  
                            "West Dutch", "East Dutch", "Daytona", "Sandy Cove"))

## format community matrix w/ specific sites: remove Sandy Cove
dat <- filter(dat, SITE %in% c("NavFac", "WestEnd Urchin", "WestEnd Kelp",  
                               "West Dutch", "East Dutch", "Daytona"))

## clear Sandy Cove info if only analyzing 6 main sites
dat<-droplevels(dat)

## remove metadata, creating community matrix of 14 spp
dat_Comm <- dat[,8:21]

## log transform 
log_Comm <- log10(dat_Comm[,1:14]+1)
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Non-Metric Multidimensional Scaling
## run NMDS 
ord <- metaMDS(comm = dat_Comm, distance="bray", k=2, min = 10, trymax=250, autotransform = F, 
               #smin = .180, sratmax = 0.9999999,
               wascores = TRUE)

## save full ordination 
save(ord, file = "ord_22March2020_noLOG.rda")

## visualize stress, check plot, xy coordinates 
stressplot(ord)
plot(ord)
scores(ord)

## coodinates saved as data frame and saved to project 
NMDS_coords<-as.data.frame(scores(ord))
#save(NMDS_coords, file="NMDS_coords")

## bind NMDS x, y, to community data 
newdf <- cbind(NMDS_coords,dat)
newdf<-droplevels(newdf)

## overlay correlation with log species  
dist <- ord$dist
ord.points <- postMDS(ord$points, dist)
scores <- wascores(ord.points, dat_Comm)     #Fall or no_SC
scores

## add species labels
plot(ord)
points(ord, display = "species")
text(ord, display = "species", col = "red")

## correlation with environmentalal variables  
ord.fit <- envfit(ord, dat$RELIEF)      #Fall or no_SC
plot(ord.fit)
ord.fit
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#### Calculate distances between points in species space for velocity calculation 
## order data (and later repeated to confirm calculations proceeded properly)
newdf <- newdf[order(newdf$PERIOD),]
newdf <- newdf[order(newdf$SITE),]
newdf <- newdf[order(newdf$TRANSECT),]

## create new DIST variable, calculate distance 
n <- 1973
newdf$DIST<-NA
newdf$DIST[2:n] <- sqrt((newdf$NMDS1[2:n] - newdf$NMDS1[1:n-1]) ^ 2 + 
                          (newdf$NMDS2[2:n] - newdf$NMDS2[1:n-1]) ^ 2)

## set initial time point to NA as we cannot calculate distance  
newdf$DIST[newdf$PERIOD == "1"] <- "NA"
newdf$DIST <- as.numeric(newdf$DIST)

## calculate distance traveled exclusively on X-axis, i.e., 1-dimensional community shift
newdf$X_dist[2:n] <- (newdf$NMDS1[2:n] - newdf$NMDS1[1:n])
newdf$X_dist[newdf$PERIOD == "1"] <- "NA"
newdf$X_dist <- as.numeric(newdf$X_dist)
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Calculate days in-between sequential sample points for velocity calculation 
## convert date into R friendly format 
newdf$date <- as.Date(newdf$date, "%m/%d/%Y")

## setup and population w/ NA column for date differences 
newdf$date_diff <- NA

## group data by ID (unique identifier for each transect, i.e., 30 unique IDs), and calculate days elapsed between individual sample points 
newdf <- newdf %>% 
  group_by(ID) %>%
  #arrange(date) %>% 
  mutate(date_diff = date - lag(date, default = first(date)))


## first convert days into character, and then into numeric 
newdf$days <- as.character(newdf$date_diff)
newdf$days <- as.numeric(newdf$days)
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
### Calculate Velocities of community shift and requisite spatial information 
## calculate Velocities by dividing 2-dimensional distances in ordination space by # of days between sampling 
newdf$Vel <- (newdf$DIST/newdf$days)

## calculate velocity in 1-dimensional space
newdf$X_Vel <- (newdf$X_dist/newdf$days)

## calculate midpoint on Axis-1 (x-axis) in order to plot velocities 
newdf$MidPt[2:n] <- ((newdf$NMDS1[1:n] + newdf$NMDS1[2:n])/2)

## set all 1st period sample points to NA, as there is no accompanying velocity to plot 
newdf$MidPt[newdf$PERIOD == "1"] <- "NA"

## drop all rows with NA (the first observation from a transect that doesn't have a distance or velocity observation)
newdf <- newdf[complete.cases(newdf),]

## save .csv with spreadsheet
write.csv(newdf,'SNI_subtidal_swath_NMDS_coordinates_noLOG.csv')
#### ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








  
