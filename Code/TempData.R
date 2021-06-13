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
library(magrittr)

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("SiteSpecificTemps.csv", header = TRUE)


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
                 legend.position = c(0.825, 0.85))

graphics.off()
windows(h=5,w=8, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data formattting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reorder sites for plotting
dat$Site <- factor(dat$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))
dat$logTemp <- log(dat$DegC)                   

## subset by site
NF <- filter(dat, Site %in% c("NavFac"))
WE <- filter(dat, Site %in% c("WestEnd"))
Day <- filter(dat, Site %in% c("Daytona"))
ED <- filter(dat, Site %in% c("EastDutch"))
ED <- na.omit(ED)

## add SU to ease plotting the xaxis 
NF$SU <- seq(from=1, to=length(NF$DegC), by=1)
WE$SU <- seq(from=1, to=length(WE$DegC), by=1)
Day$SU <- seq(from=1, to=length(Day$DegC), by=1)
ED$SU <- seq(from=1, to=length(ED$DegC), by=1)

newDat <- rbind(NF, WE, Day, ED)


## Plot entire time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color pal matching ms site colors. NF, blend(WEU & WEK), Day, ED
pal_sites <- c("#BE2625", "#E88600", "#608341", "#00536f") 
#pal_sites2 <- c("#BE2625", "#E88600", "#6f006b", "#00536f") #BE2625 used for tint and shade creation

graphics.off()
windows(w=12,h=4,record=TRUE)


## all sites 
p1 <- ggplot(newDat, aes(SU, DegC, color=Site, group=Site)) + my.theme +
  geom_line(alpha=.5) + ylab("Degrees Celcius") + 
  scale_colour_manual(values=pal_sites) +
  scale_x_continuous(n.breaks = 17, labels=c("remove", "Nov 2015", 
                                             "Mar 2016", "June 2016", "Oct 2016", "Jan 2017",
                                             "May 2017", "Aug 2017", "Dec 2017", "Mar 2018",
                                             "July 2018", "Oct 2018", "Feb 2019", "May 2019",
                                             "Sep 2019", "")) +
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        legend.position = c(0.825, 0.85))

print(p1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## visualize kernal densities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
graphics.off()
windows(w=8,h=5,record=TRUE)

## ind sites
p2 <- ggplot(dat, aes(DegC, fill=Site)) +
  geom_density(binwidth = 0.15, color="black") +  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  facet_wrap(~Site)
print(p2)


## overlapping kernal densities
p3 <- ggplot(dat, aes(DegC, fill=Site)) +
  geom_density(position="identity", binwidth = 0.15, color="black", alpha=0.3) +  my.theme +
  scale_fill_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") +
  guides(fill=guide_legend(order=1))
print(p3)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS two-sample tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
## END acf and ccf ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## inverse ecdf plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## calcuate empirical cumulative density function and extract and sort values
NF_ecdf <- data.frame(x=unique(NF$DegC), 
                      y=ecdf(NF$DegC)(unique(NF$DegC))*length(NF$DegC))
WE_ecdf <- data.frame(x=unique(WE$DegC), 
                      y=ecdf(WE$DegC)(unique(WE$DegC))*length(WE$DegC))
ED_ecdf <- data.frame(x=unique(ED$DegC), 
                      y=ecdf(ED$DegC)(unique(ED$DegC))*length(ED$DegC))
Day_ecdf <- data.frame(x=unique(Day$DegC), 
                       y=ecdf(Day$DegC)(unique(Day$DegC))*length(Day$DegC))


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
pal_sites <- c("#BE2625", "#E88600", "#608341", "#00536f") #BE2625 used for tint and shade creation


## plot all 
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=Site)) + my.theme +
  geom_line(lwd=1, alpha=.8) +
  scale_color_manual(values=pal_sites) +
  xlab("Temperature, degrees Celcius") + ylab("empirical probability") 
  #ggtitle("probability of water temperature being equal to or warmer than the corresponding (x-axis) temp")

print(p1)
## END inverse ecdf plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


