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


## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data formattting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reorder sites for plotting
dat$Site <- factor(dat$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))
                   

## subset by site
NF <- filter(dat, Site %in% c("NavFac"))
WE <- filter(dat, Site %in% c("WestEnd"))
Day <- filter(dat, Site %in% c("Daytona"))
ED <- filter(dat, Site %in% c("EastDutch"))





## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color pal matching ms site colors. NF, blend(WEU & WEK), Day, ED
pal_sites <- c("#BE2625", "#E88600", "#608341", "#00536f") #BE2625 used for tint and shade creation


## graphing window
graphics.off()
windows(w=12,h=6,record=TRUE)


## one site at a time 
p1 <- ggplot(Day, aes(DateTime, DegC)) +
  geom_line()
print(p1)


## all sites 
p5 <- ggplot(dat, aes(DateTime, DegC, color=Site)) +
  geom_line() +
  scale_colour_manual(values=pal_sites) +
  facet_wrap(~ Site)
print(p5)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## visualize frequency histograms 
p6 <- ggplot(dat, aes(DegC, fill=Site)) +
  geom_histogram(binwidth = 0.15, color="black") + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  ggtitle("site specific frequency histograms for 2016-2019 temperature data") +
  facet_wrap(~Site)
print(p6)



## transparent and overlapping histograms.
p6 <- ggplot(dat, aes(DegC, fill=Site)) +
  geom_histogram(position="identity", binwidth = 0.15, color="black", alpha=0.75) + 
  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  ggtitle("overlapping histograms for 2016-2019 temperature data") +
  guides(fill=guide_legend(order=1))
print(p6)



## Degrees C as a numeric vector
NF_num <- NF$DegC
WE_num <- WE$DegC
Day_num <- Day$DegC
ED_num <- ED$DegC
ED_num <- na.omit(ED_num)



## perform KS tests 
ks.test(NF_num, ED_num)
ks.test(NF_num, WE_num)
ks.test(NF_num, Day_num)
ks.test(WE_num, ED_num)
ks.test(WE_num, Day_num)
ks.test(Day_num, ED_num)


## simulate minimally to get a sense of the behavior of D
x <- rnorm(3000, mean = 0, sd = 2)
y <- rnorm(3000, mean = 10, sd = 2)
xy <- ks.test(x,y)
## END KS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## wavelet analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## data frames for wavelet analysis
NF_temp <- as.data.frame(NF$DegC)
WE_temp <- as.data.frame(WE$DegC)
ED_temp <- as.data.frame(ED$DegC)
Day_temp <- as.data.frame(Day$DegC)


## there are 9 NAs in East Dutch (sensor issues?) remove them for wavelet 
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




## visualize wavelet analysis 
graphics.off()
windows(w=12,h=8,record=TRUE)


## unsure why par is not working here? 
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
## END wavelet analysis ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## auto-correlation and cross-correlation through time ~~~~~~~~~~~~~~~~~~~~~~~~~
auto_cor <- acf(NF_temp, lag.max = 16000, 
                main="auto-correlation within NavFac temp data", 
                xlab="lag in units of 1 hr; purple line = 6 month lag, red line = 1 year lag",
                ylab="auto-correlation")
abline(v=8760, col="red")
abline(v=4380, col="darkorchid")


## cross-correlation through time between two time series 
cross_cor <- ccf(NF_temp, ED_temp, type = "correlation", 
                 lag.max = 16000,
                 main="cross-correlation through time between NacFac and EastDutch temp data",
                 xlab="lag in units of 1 hr; purple line = 6 month lag, red line = 1 year lag",
                 ylab="cross-correlation")
abline(v=8760, col="red")
abline(v=4380, col="darkorchid")



