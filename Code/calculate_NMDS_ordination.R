## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## non-metric multidimensional scaling (NMDS) analyses of SNI subtidal data 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())


## load packages (if packages not installed: install.packages("librarynamehere"))
library(vegan)
library(MASS)
library(fitdistrplus)
library(tidyverse)
library(plyr)


## set path 
setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")


## data upon which to calculate Bray-Curtis and to perform NMDS analysis 
dat <- read.csv("SNI_SubtidalSwath_OriginalDataFile.csv", header = TRUE)


## load ordination (results of a previously run ordination, so one does not have to take large amounts of time to rerun analysis) 
load("ord_22March2020.rda")


## rugosity data 
rugosity <- read.csv("SubstrateRugosity.csv", header = TRUE)
## END startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## process rugosity values and add to community matrix ~~~~~~~~~~~~~~~~~~~~~~~~~
## calculate mean, sd, se of rugosity
rugosity.summary <- ddply(rugosity, c("SITE", "TRANSECT", "ID"), summarise, 
                          N = length(RELIEF),
                          mean_rug = mean(RELIEF),
                          sd_rug = sd(RELIEF),
                          se_rug = sd_rug / sqrt(N))


## filter out sandy cove
rugosity.summary <- filter(rugosity.summary, SITE %in% c("NavFac","West End Kelp","West End Urchin",
                                                        "Daytona","East Dutch","West Dutch"))


## join rugosity values to ordination data 
dat <- dplyr::inner_join(dat, rugosity.summary, by="ID")


## fix names of original file columns (altered from table join)
names(dat)[4]<-"SITE"
names(dat)[6]<-"TRANSECT"


## delete redunant columns (added from table join)
dat <- dat[, !(colnames(dat) %in% c("SITE.y","TRANSECT.y","N"))]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## set site as factor and rename 
dat$SITE <- factor(dat$SITE, levels=c("1", "2", "3", "4",
                                      "5", "6"), 
                   labels=c("NavFac", "WestEnd Urchin", "WestEnd Kelp",  
                            "West Dutch", "East Dutch", "Daytona"))


## format community matrix w/ specific sites: remove Sandy Cove
dat <- filter(dat, SITE %in% c("NavFac", "WestEnd Urchin", "WestEnd Kelp",  
                               "West Dutch", "East Dutch", "Daytona"))


## remove metadata, creating community matrix of 14 spp
dat_Comm <- dat[,7:20]

## log transform 
log_Comm <- log10(dat_Comm[,1:14]+1)
## END data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Non-Metric Multidimensional Scaling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#ord <- metaMDS(comm = log_Comm, distance="bray", k=2, autotransform = F, wascores = TRUE)


## save full ordination 
#save(ord, file = "ord_date.rda")


## visualize stress, check plot, xy coordinates 
#stressplot(ord)
#plot(ord)
#scores(ord)


## overlay correlation with log species  
#dist <- ord$dist
#ord.points <- postMDS(ord$points, dist)
#scores <- wascores(ord.points, dat_Comm)     
#scores


## add species labels
#plot(ord)
#points(ord, display = "species")
#text(ord, display = "species", col = "red")


## correlation with environmentalal variables  
#ord.fit <- envfit(ord, dat$mean_rug)      
#plot(ord.fit)
#ord.fit


## coodinates saved as data frame and bound to community / site dataframe
bind.coordinates <- function(ordination, data){
  t1 <- as.data.frame(scores(ordination))
  t2 <- cbind(t1, data)
  return(t2)
  }


## invoke function
dat <- bind.coordinates(ord, dat)
## END NMDS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## Calculate distances between points in species space for velocity calculation  
## order data (and later repeated to confirm calculations proceeded properly)
dat <- dat[order(dat$PERIOD),]
dat <- dat[order(dat$SITE),]
dat <- dat[order(dat$TRANSECT),]

## create new DIST variable, calculate distance 
n <- 1973
dat$DIST<-NA
dat$DIST[2:n] <- sqrt((dat$NMDS1[2:n] - dat$NMDS1[1:n-1]) ^ 2 + 
                          (dat$NMDS2[2:n] - dat$NMDS2[1:n-1]) ^ 2)

## set initial time point to NA as we cannot calculate distance  
dat$DIST[dat$PERIOD == "1"] <- "NA"
dat$DIST <- as.numeric(dat$DIST)

## calculate distance traveled exclusively on X-axis, i.e., 1-dimensional community shift
dat$X_dist[2:n] <- (dat$NMDS1[2:n] - dat$NMDS1[1:n])
dat$X_dist[dat$PERIOD == "1"] <- "NA"
dat$X_dist <- as.numeric(dat$X_dist)
## END distance calculation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Calculate days in-between sequential sample points for velocity calculation ~
## convert date into R friendly format 
dat$date <- as.Date(dat$date, "%m/%d/%Y")


## setup and population w/ NA column for date differences 
dat$date_diff <- NA


## group data by ID (unique identifier for each transect, i.e., 30 unique IDs), and calculate days elapsed between individual sample points 
dat <- dat %>% 
  group_by(ID) %>%
  #arrange(date) %>% 
  mutate(date_diff = date - lag(date, default = first(date)))


## first convert days into character, and then into numeric 
dat$days <- as.character(dat$date_diff)
dat$days <- as.numeric(dat$days)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Calculate Velocities of community shift and requisite spatial information ~~~
## calculate Velocities by dividing 2-dimensional distances in ordination space by # of days between sampling 
dat$Vel <- (dat$DIST/dat$days)


## calculate velocity in 1-dimensional space
dat$X_Vel <- (dat$X_dist/dat$days)


## calculate midpoint on Axis-1 (x-axis) in order to plot velocities 
dat$MidPt[2:n] <- ((dat$NMDS1[1:n] + dat$NMDS1[2:n])/2)


## set all 1st period sample points to NA, as there is no accompanying velocity to plot 
dat$MidPt[dat$PERIOD == "1"] <- "NA"


## drop all rows with NA (the first observation from a transect that doesn't have a distance or velocity observation)
dat <- dat[complete.cases(dat),]


## set all rows with < 0 days to 0 (these are the Period = 1 sample points)
dat$days[dat$days <= "0"] <- "0"


## remove these Period = 1 transects 
dat = filter(dat, days > "0", )


## save .csv with spreadsheet
write.csv(dat,'SNI_subtidal_swath_NMDS_coordinates_2Removed.csv')
## END velocity calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
