## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## SNI temperature data Fig S5 in SubstrateComplexity ~~~~~~~~~~~~~~~~~~~~~~~ ##
## Modified May 25th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
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
filter.f <- function(x){
  t1 <- filter(dat, Site %in% c(x))
  SU <- seq(from=1, to=length(t1$DegC), by=1)
  out <- return(as.data.frame(cbind(SU, t1)))
}

NF <- filter.f("NavFac")
WE <- filter.f("WestEnd")
Day <- filter.f("Daytona")
ED <- na.omit(filter.f("EastDutch"))
newDat <- rbind(NF, WE, Day, ED)
## END data formatting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot entire time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color pal matching ms site colors. NF, blend(WEU & WEK), Day, ED
cols <- scale_color_manual(values=c("#BE2625", "#E88600", "#608341", "#00536f")) 

x.breaks <- scale_x_continuous(n.breaks = 17, 
                               labels=c("remove", "Nov 2015", "Mar 2016", "June 2016", "Oct 2016", "Jan 2017",
                                         "May 2017", "Aug 2017", "Dec 2017", "Mar 2018", "July 2018", "Oct 2018", 
                                         "Feb 2019", "May 2019","Sep 2019", ""))
  

## all sites 
p1 <- ggplot(newDat, aes(SU, DegC, color=Site, group=Site)) +  
  geom_line(alpha=.5) + ylab("Degrees Celcius") + my.theme + cols + x.breaks + 
  theme(axis.text.x = element_text(angle=45, hjust=1), axis.title.x = element_blank(),
        legend.position = c(0.825, 0.85))

graphics.off()
windows(w=12,h=4,record=TRUE)
print(p1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## visualize kernal densities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
## standardize sample size among the four sites 
len.f <- function(x1, x2, x3, x4){
l1 <- length(x1[,4])
l2 <- length(x2[,4])
l3 <- length(x3[,4])
l4 <- length(x4[,4])
sample.size <- min(l1, l2, l3, l4)
out <- return(sample.size)
}

## calculate mimimum sample size used to rarify data 
sample.size <- len.f(NF, WE, Day, ED)


## rarify data and return data data frame
sample.f <- function(x, site){
  out <- as.data.frame(sample(x[,3], sample.size, replace=F))
  names(out)[1] <- "DegC"
  out$site <- site
  return(out)
  }

NavFac <- sample.f(NF, "NavFac")
WestEnd <- sample.f(WE, "WestEnd")
Daytona <- sample.f(Day, "Daytona")
EastDutch <- sample.f(ED, "EastDutch")

newDat <- rbind(NavFac, WestEnd, Daytona, EastDutch)


## overlapping kernal densities
p3 <- ggplot(dat, aes(DegC, fill=Site)) +
  geom_density(position="identity", color="black", alpha=0.3) +  my.theme + cols + 
  xlab("Temperature, degrees Celcius") + guides(fill=guide_legend(order=1))

graphics.off()
windows(w=8,h=5,record=TRUE)
print(p3)
## END kernal density visualization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## KS two-sample tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Degrees C as a numeric vector
NF_num <- NavFac$DegC
WE_num <- WestEnd$DegC
Day_num <- Daytona$DegC
ED_num <- EastDutch$DegC


## perform KS tests 
ks.test(NF_num, ED_num)
ks.test(NF_num, WE_num)
ks.test(NF_num, Day_num)
ks.test(WE_num, ED_num)
ks.test(WE_num, Day_num)
ks.test(Day_num, ED_num)
## END KS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## inverse ecdf plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## t1 = calculate empirical cumulative density function values
## t2 = rescale extracted cdf values to 0-1 scale
## t3 = take the inverse of a cdf, such that the p(x) > or = log wave event
ecdf.f <- function(dat, site.name){
  t1 <- data.frame(x=unique(dat$DegC), y=ecdf(dat$DegC)(unique(dat$DegC))*length(dat$DegC))
  t2 <- scale(t1$y, center=min(t1$y), scale=diff(range(t1$y)))
  t3 <- ((t2 - max(t2)) * (-1))
  t4 <- cbind(t1, t2, t3)
  names(t4)[1:4]<-c("x", "unscaled_y", "y", "inv_y")
  t4$Site <- site.name 
  out <- return(t4)
}


## perform calculations
NF_ecdf <- ecdf.f(NF, "NavFac")
WE_ecdf <- ecdf.f(WE, "WestEnd")
Day_ecdf <- ecdf.f(Day, "Daytona")
ED_ecdf <- ecdf.f(ED, "EastDutch")


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat_ecdf <- rbind(NF_ecdf, WE_ecdf, ED_ecdf, Day_ecdf)
dat_ecdf$Site <- factor(dat_ecdf$Site, levels=c("NavFac", "WestEnd", "Daytona", "EastDutch"))


## plot all 
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=Site)) + my.theme + cols +
  geom_line(lwd=1, alpha=.8) + xlab("Temperature, degrees Celcius") + ylab("empirical probability") 
print(p1)
## END inverse ecdf plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
