## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Visualization / Analyses of 2017 Quadrat Survey around San Nicolas Island 
## zhr -- May 23rd 2021 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## clear Global Environment 
rm(list = ls())

## load packages (if packages not installed: install.packages("librarynamehere"))
library(ggplot2)
library(dplyr)
library(vegan)
library(MASS)
library(fitdistrplus)
library(pryr)
library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
library(beeswarm)

## close out window tab & open custom window
graphics.off()
windows(h=6,w=9, record=TRUE)

## set working directory
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")

## load data frame
dat <- read.csv("QuadSurvey2017.csv", header = TRUE)

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data organization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## only include new Zones and not quadrats sampled along the index site
dat <- filter(dat, Zone %in% c("1","2","3")) #,"5","6"))


## Create new column with sites broken down by relief to expedite plotting 
dat$Ind <- NA
dat$Ind[dat$Relief=="Low" & dat$Site=="NavFac"]<-"NavFac_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="NavFac"]<-"NavFac_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="East Dutch"]<-"EastDutch_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="East Dutch"]<-"EastDutch_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="West Dutch"]<-"WestDutch_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="West Dutch"]<-"WestDutch_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="West End"]<-"WestEnd_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="West End"]<-"WestEnd_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="Daytona"]<-"Daytona_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="Daytona"]<-"Daytona_High"




## reorder data (for plotting)
dat$Ind <- factor(dat$Ind, levels = c("NavFac_Low", "NavFac_High", 
                                        "WestEnd_Low", "WestEnd_High",
                                        "Daytona_Low", "Daytona_High",
                                        "EastDutch_Low", "EastDutch_High",
                                        "WestDutch_Low", "WestDutch_High"))



## custom colors
color_pal <- c("#9D1309","#0E808B")



## split up data types; transform coverage to 0-1 scale
metaData <- dat[,c(1:9,31)]
percent <- dat[,10:20] / 100
abundance <- dat[,21:30]
dat <- cbind(metaData,percent,abundance)
## End data organization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## visualize verticality of quadrats
beeswarm(Verticality ~ Ind, data = dat, main = 'Quadrat elevation above seafloor',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Vertical height (m) above seafloor", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)





## close out window tab & open custom window
graphics.off()
windows(h=10,w=15, record=TRUE)
par(mfrow=c(2,2))

## Select percent coverage results ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## cucumbers 
beeswarm(Cucumber ~ Ind, data = dat, main = '(a) Cucumber % cover',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Cucumber % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(Cucumber ~ Ind, data = dat, probs=.5, add = TRUE)


## fleshy red algae
beeswarm(FleshyRed ~ Ind, data = dat, main = '(b) Fleshy red algae % cover',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="FleshyRed % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(FleshyRed ~ Ind, data = dat, probs=.5, add = TRUE)


## Articulated coralline algae
beeswarm(ArtCor ~ Ind, data = dat, main = '(c) Articulated coralline % cover',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Articulated Coralline % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(ArtCor ~ Ind, data = dat, probs=.5, add = TRUE)


## Encrusting coralline algae
beeswarm(EncrustCor ~ Ind, data = dat, main = '(d) Encrusting coralline % cover',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Encrusting Coralline % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(EncrustCor ~ Ind, data = dat, probs=.5, add = TRUE)
## End % cover ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Select Abundance data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## close out window tab & open custom window
graphics.off()
windows(h=10,w=15, record=TRUE)
par(mfrow=c(2,2))


beeswarm(Cystoseira ~ Ind, data = dat, main = '(a) Stephanocystis density',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Cystoseira abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topright", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(Cystoseira ~ Ind, data = dat, probs=.5, add = TRUE)




beeswarm(LamSpp ~ Ind, data = dat, main = '(b) Laminaria spp. density',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Laminaria spp. abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(LamSpp ~ Ind, data = dat, probs=.5, add = TRUE)




beeswarm(PurpUrch ~ Ind, data = dat, main = '(c) Purple urchin density',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Purple Urchin abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topright", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(PurpUrch ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(RedUrch ~ Ind, data = dat, main = '(d) Red urchin density',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Red Urchin abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(RedUrch ~ Ind, data = dat, probs=.5, add = TRUE)
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##









