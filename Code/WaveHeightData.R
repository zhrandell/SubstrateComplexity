## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## SNI wave height data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Modified June 1st 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(ggplot2)
library(dplyr)
library(vegan)
library(MASS)
library(fitdistrplus)
library(gridExtra)
library(RColorBrewer)
library(ggpubr)
library(egg)
library(gtable)
library(grid)
library(magick)
library(magrittr)
library(here)
library(mclust, quietly=TRUE)
library(pryr)
library(ggbeeswarm)
library(WaveletComp)
library(devtools)


setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("WaveHeights.csv", header = TRUE)


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





## data formattting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reorder sites for plotting
dat$Site <- factor(dat$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))
                   

## restrict data to sample points where all sensors are simultaneously deployed
dat <- filter(dat, SU > 3364864 & SU < 28870656)


## log transform
dat$logWV <- log(dat$WaveHeight)


## subset by site
NF <- filter(dat, Site %in% c("NavFac"))
WE <- filter(dat, Site %in% c("WestEnd"))
Day <- filter(dat, Site %in% c("Daytona"))
ED <- filter(dat, Site %in% c("EastDutch"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color pal matching ms site colors. NF, blend(WEU & WEK), Day, ED
pal_sites <- c("#BE2625", "#E88600", "#608341", "#00536f") #BE2625 used for tint and shade creation
#pal_sites2 <- c("#BE2625", "#E88600", "#6f006b", "#00536f") #BE2625 used for tint and shade creation

#6f006b
## graphing window
graphics.off()
windows(w=12,h=4,record=TRUE)


## all sites 
p5 <- ggplot(dat, aes(SU, WaveHeight, color=Site)) +
  geom_line(alpha=.7) +
  scale_colour_manual(values=pal_sites) +
  scale_x_continuous(n.breaks=12, labels=c("null", 
                                           "Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug", "Sep",
                                           "")) +
  my.theme + theme(axis.title.x=element_blank()) + ylab("Wave height (meters)")
print(p5)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## visualize frequency histograms 
p6 <- ggplot(dat, aes(WaveHeight, fill=Site)) +
  geom_histogram(binwidth = 0.015) + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Wave height, meters") +
  ggtitle("site specific frequency histograms for one year of wave height data") +
  facet_wrap(~Site)
print(p6)


p6 <- ggplot(dat, aes(WaveHeight, fill=Site)) +
  geom_density(binwidth = 0.15, color="black") + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  ggtitle("site specific frequency histograms for 2016-2019 temperature data") +
  facet_wrap(~Site)
print(p6)




## transparent and overlapping histograms.
p6 <- ggplot(dat, aes(WaveHeight, fill=Site)) +
  geom_histogram(position="identity", binwidth = 0.015, alpha=0.4) + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  ggtitle("overlapping histograms for 2016-2019 temperature data") +
  guides(fill=guide_legend(order=1))
print(p6)


## graphing window
graphics.off()
windows(w=8,h=5,record=TRUE)

## use logWV for logged data, and WaveHeight for "normal" scale (meters)
p7 <- ggplot(dat, aes(logWV, fill=Site)) +
  geom_density(position="identity", color="black", alpha=0.3) + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Wave height (log meters)") +
  #ggtitle("overlapping kernal densities for one year of wave height data") +
  guides(fill=guide_legend(order=1))
print(p7)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## does wave height data follow a power relationship? ~~~~~~~~~~~~~~~~~~~~~~~~~~
## calcuate empirical cumulative density function and extract and sort values
NF_ecdf <- data.frame(x=unique(NF$logWV), 
                     y=ecdf(NF$logWV)(unique(NF$logWV))*length(NF$logWV))
WE_ecdf <- data.frame(x=unique(WE$logWV), 
                     y=ecdf(WE$logWV)(unique(WE$logWV))*length(WE$logWV))
ED_ecdf <- data.frame(x=unique(ED$logWV), 
                     y=ecdf(ED$logWV)(unique(ED$logWV))*length(ED$logWV))
Day_ecdf <- data.frame(x=unique(Day$logWV), 
                     y=ecdf(Day$logWV)(unique(Day$logWV))*length(Day$logWV))


## rescale extracted cdf values to 0-1 scale
NF_ecdf$y <- scale(NF_ecdf$y, center=min(NF_ecdf$y), scale=diff(range(NF_ecdf$y)))
WE_ecdf$y <- scale(WE_ecdf$y, center=min(WE_ecdf$y), scale=diff(range(WE_ecdf$y)))
ED_ecdf$y <- scale(ED_ecdf$y, center=min(ED_ecdf$y), scale=diff(range(ED_ecdf$y)))
Day_ecdf$y <- scale(Day_ecdf$y, center=min(Day_ecdf$y), scale=diff(range(Day_ecdf$y)))


## take the inverse of a cdf, such that the p(x) > or = log wave event
NF_ecdf$inv_y <- ((NF_ecdf$y - max(NF_ecdf$y)) * (-1))
WE_ecdf$inv_y <- ((WE_ecdf$y - max(WE_ecdf$y)) * (-1))
ED_ecdf$inv_y <- ((ED_ecdf$y - max(ED_ecdf$y)) * (-1))
Day_ecdf$inv_y <- ((Day_ecdf$y - max(Day_ecdf$y)) * (-1))


## add site name to new data frames
NF_ecdf$Site <- "NavFac"
WE_ecdf$Site <- "WestEnd"
ED_ecdf$Site <- "EastDutch"
Day_ecdf$Site <- "Daytona"


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat_ecdf <- rbind(NF_ecdf, WE_ecdf, ED_ecdf, Day_ecdf)
dat_ecdf$Site <- factor(dat_ecdf$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))


## plot all 
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=Site)) + my.theme +
  geom_line(lwd=1, alpha=.8) +
  scale_color_manual(values=pal_sites) +
  xlab("Wave height (log meters)") + ylab("inverse empirical CDF")  
  #ggtitle("probability of a wave event of equal or greater size relative to log wave height")

print(p1)
## End script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~








