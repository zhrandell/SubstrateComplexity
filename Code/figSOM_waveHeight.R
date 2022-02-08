## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## SNI wave height data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Modified June 1st 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
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
## END initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data formattting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reorder sites for plotting
dat$Site <- factor(dat$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))
                   

## restrict data to sample points where all sensors are simultaneously deployed
dat <- filter(dat, SU > 3364864 & SU < 28870656)


## log transform
dat$logWV <- log10(dat$WaveHeight)


## filter by site
filter.site <- function(x){filter(dat, Site %in% c(x))}
NF <- filter.site("NavFac")
WE <- filter.site("WestEnd")
Day <- filter.site("Daytona")
ED <- filter.site("EastDutch")
## END data formatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cols <- scale_color_manual(values=c("#BE2625", "#E88600", "#608341", "#00536f")) 
x.lab.ts <- scale_x_continuous(n.breaks=12,labels=c("null","Nov","Dec","Jan","Feb","Mar","Apr","May","June","July","Aug", "Sep","")) 
x.lab.kd <- scale_x_continuous(labels=c("null", "0.40", "1", "2.5", "6.3")) 
## END custom params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot time series and kernal densities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot time series of untransformed wave height
p1 <- ggplot(dat, aes(SU, WaveHeight, color=Site)) + geom_line(alpha=.7) + cols + x.lab.ts + 
  my.theme + theme(axis.title.x=element_blank()) + ylab("Wave height (meters)")


## plot kernal density use log base-10 transformed wave data
p2 <- ggplot(dat, aes(logWV, fill=Site)) + geom_density(position="identity", color="black", alpha=0.3) + 
  my.theme + x.lab.kd + cols + xlab("Wave height (meters)") + guides(fill=guide_legend(order=1))


## plot wave height time series
graphics.off()
windows(w=12,h=4,record=TRUE)
print(p1)


## plot wave height kernal densities
windows(w=8,h=5,record=TRUE)
print(p2)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS tests to test kernal density distributions between-sites ~~~~~~~~~~~~~~~~~
NF_num <- NF$logWV
WE_num <- WE$logWV
Day_num <- Day$logWV
ED_num <- ED$logWV

## perform KS tests 
ks.test(NF_num, ED_num)
ks.test(NF_num, WE_num)
ks.test(NF_num, Day_num)
ks.test(WE_num, ED_num)
ks.test(WE_num, Day_num)
ks.test(Day_num, ED_num)
## END KS tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## calcuate empirical cumulative density function and extract/sort/process values 
ecdf.f <- function(dat, site.name){
  t1 <- data.frame(x=unique(dat$logWV), y=ecdf(dat$logWV)(unique(dat$logWV))*length(dat$logWV))
  t2 <- scale(t1$y, center=min(t1$y), scale=diff(range(t1$y)))
  t3 <- ((t2 - max(t2)) * (-1))
  t4 <- cbind(t1, t2, t3)
  names(t4)[1:4] <- c("x", "unscaled_y", "y", "inv_y")
  t4$Site <- site.name
  out <- return(t4)
}


## apply ecdf.f() function to data
NF.e <- ecdf.f(NF, "NavFac")
WE.e <- ecdf.f(WE, "WestEnd")
ED.e <- ecdf.f(ED, "EastDutch")
Day.e <- ecdf.f(Day, "Daytona")


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat.e <- rbind(NF.e, WE.e, ED.e, Day.e)
dat.e$Site <- factor(dat.e$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))


## plot all 
p3 <- ggplot(dat.e, aes(x, inv_y, color=Site)) + my.theme + geom_line(lwd=1, alpha=.8) +
  xlab("Wave height (meters)") + ylab("inverse empirical CDF") + cols + x.lab.kd
print(p3)
## END ecdf plot~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
