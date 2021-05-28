## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## SNI temperature data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Modified May 25th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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


setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("SiteSpecificTemps.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## subset by site
NF <- filter(dat, Site %in% c("NavFac"))
WE <- filter(dat, Site %in% c("WestEnd"))
Day <- filter(dat, Site %in% c("Daytona"))
ED <- filter(dat, Site %in% c("EastDutch"))



## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(w=24,h=12,record=TRUE)


## one site at a time 
p1 <- ggplot(Day, aes(DateTime, DegC)) +
  geom_line()
print(p1)


## all sites 
p5 <- ggplot(dat, aes(DateTime, DegC, color=Site)) +
  geom_line() +
  facet_wrap(~ Site)
print(p5)





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## data frames for wavelet analysis
NF_temp <- as.data.frame(NF$DegC)
WE_temp <- as.data.frame(WE$DegC)
ED_temp <- as.data.frame(ED$DegC)
Day_temp <- as.data.frame(Day$DegC)


ED_temp <- na.omit(ED_temp)



## wavelet analysis 

NF_wavelet = analyze.wavelet(NF_temp, loess.span = 0,
                             dt = 1, dj = 1/20,
                             lowerPeriod = 1,
                             upperPeriod = 16000,
                             make.pval = F, n.sim = 10)


WE_wavelet = analyze.wavelet(WE_temp, loess.span = 0,
                             dt = 1, dj = 1/20,
                             lowerPeriod = 1,
                             upperPeriod = 16000,
                             make.pval = F, n.sim = 10)


ED_wavelet = analyze.wavelet(ED_temp, loess.span = 0,
                              dt = 1, dj = 1/20,
                              lowerPeriod = 1,
                              upperPeriod = 16000,
                              make.pval = F, n.sim = 10)


Day_wavelet = analyze.wavelet(Day_temp, loess.span = 0,
                        dt = 1, dj = 1/20,
                        lowerPeriod = 1,
                        upperPeriod = 16000,
                        make.pval = F, n.sim = 10)





graphics.off()
windows(w=24,h=12,record=TRUE)

par(mfrow=c(2,2))


NavFac <- wt.image(NF_wavelet, main = "NavFac", color.key = "quantile", n.levels = 250,
         legend.params = list(lab = "wavelet power levels"))


WestEnd <- wt.image(WE_wavelet, main = "WestEnd", color.key = "quantile", n.levels = 250, 
                   legend.params = list(lab = "wavelet power levels"))


EastDutch <- wt.image(ED_wavelet, main = "EastDutch", color.key = "quantile", n.levels = 250, 
                    legend.params = list(lab = "wavelet power levels"))


Daytona <- wt.image(Day_wavelet, main = "Daytona", color.key = "quantile", n.levels = 250, 
                      legend.params = list(lab = "wavelet power levels"))


dev.off()


