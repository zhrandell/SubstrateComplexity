## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Script to produce Figure 3, Urchin Densities for Substrate Complexity ms ~ ##
## updated May 25th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

## set window size
graphics.off()
windows(h=3,w=3, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"), 
                   labels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                            "Daytona", "East Dutch", "West Dutch"))


# calculate total urchin densities, algae densities (unused here), and ratios
dat$totalUrchin<-dat$StrPur+dat$MesFra
dat$totalAlgae<-dat$SteOsm+dat$LamSpp+dat$MacJuv+dat$LamJuv+dat$EisArb+dat$PteCal+dat$MacPyr+1
dat$ratio <- dat$totalUrchin/dat$totalAlgae


## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))
## end data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## URCHIN DENSITY BEESWARM PLOTS ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## check kernal density
NF_kd <- ggplot(NF, aes(NMDS1)) +
  geom_density()
print(NF_kd)


## set empty pane to house final product
NF_kd_print <- ggplot(NF, aes(NMDS1)) +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + 
  theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 


### extract vectors x and y coordinates from NMDS kernal density plot
S1 <- density(NF$NMDS1)$x
S2 <- density(NF$NMDS1)$y


## calculate max y value; calculate associated x coordinate 
max_y1 <- which.max(S2)
max_x1 <- S1[max_y1]


## calculate second max y point from kernal density (modify < value accordingly)
below <- max(S2[S1<0])
belowY<-which(S2==below)
max_x2 <- S1[belowY]


## calculate minimum y point (trough) and associated x between the two maximums
trough <- min(S2[S1 < max_x1 & S1 > max_x2])
ymin<-which(S2==trough)
min_x1<-S1[ymin]


## check plot 
ggplot(NF, aes(NMDS1)) + geom_density() + 
  geom_vline(xintercept=max_x1) +
  geom_vline(xintercept=max_x2) +
  geom_vline(xintercept=min_x1)

## create new data frames from respective community states
NF_urchin <- filter(NF, NMDS1>min_x1)
NF_mixed <- filter(NF, NMDS1<min_x1)


## create beeswarm cluster from the mixed state
n1<-ggplot(NF_mixed) +
  geom_boxplot(data=NF_mixed, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + 
  ylim(0,1600) + 
  geom_beeswarm(data=NF_mixed, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,groupOnX = NULL, cex=5, stroke=.35,priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## create beeswarm cluster from the urchin-barren state
n2<-ggplot(NF_urchin) +
  geom_boxplot(data=NF_urchin, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + 
  ylim(0,1600) + 
  geom_beeswarm(data=NF_urchin, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8, groupOnX = NULL, cex=5, stroke=.35,priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n2)$data[[1]]
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## coordinates to match labels up with kernal density placement
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotation for the figure
nf <- grobTree(text_grob("NavFac", x=.4, y=.95, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.33333, y=.045, hjust=0, size = 13, color = "black"))
B <- grobTree(text_grob("B", x=.63333, y=.045, hjust=0, size = 13, color = "black"))


## final figure
NF_out <-NF_kd_print + 
  annotation_custom(ggplotGrob(n2), xmin = (high2-width), xmax = (high2+width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n1), xmin = low2-width, xmax = low2+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines")) +
 #geom_text(x=.33333,y=.1,label="M",size=4) +
 #geom_text(x=.66666,y=.1,label="B",size=4) +
  annotation_custom(nf) +
  annotation_custom(M) +
  annotation_custom(B)

print(NF_out)
## END NavFac plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## WestEnd plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom blank plot for final figure 
WEK_kd <- ggplot(WEK, aes(NMDS1)) +
  geom_density(adjust=0.8) 
print(WEK_kd)
  
WEK_kd_print <- ggplot(WEK, aes(NMDS1)) + 
  #geom_density() +
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + 
  theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 



## extract data from kernal density
p<-ggplot_build(WEK_kd)
head(p$data[[1]], 63)
kd<-as.data.frame(p$data[[1]])


## extract x and y of NMDS kernal density
S1 <- kd$x
S2 <- kd$y


## calculate maximum y and associated x value 
max_y1 <- which.max(S2)
max_x1 <- S1[max_y1]


## second max x and y value
below <- max(S2[S1< -0.75])
belowY<-which(S2==below)
max_x2 <- S1[belowY]


## calculate min between two states 
trough <- min(S2[S1 < max_x1 & S1 > max_x2])
ymin<-which(S2==trough)
min_x1<-S1[ymin]


## ID third max (for algae state)
below2 <- max(S2[S1 > .5 & S1 < 1])
belowY2<-which(S2==below2)
max_x3 <- S1[belowY2]


## calculate 2nd min 
trough <- min(S2[S1 > max_x1 & S1 < max_x3])
ymin2<-which(S2==trough)
min_x2<-S1[ymin2]


## check plot 
WEK_kd +
  geom_vline(xintercept=max_x1) +
  geom_vline(xintercept=max_x2) +
  geom_vline(xintercept=min_x1) +
  geom_vline(xintercept=max_x3) +
  geom_vline(xintercept=min_x2) 


## new data frames
WEK_urchin <- filter(WEK, NMDS1>min_x2)
WEK_mixed <- filter(WEK, NMDS1<min_x2&NMDS1>min_x1)
WEK_algae <- filter(WEK, NMDS1<min_x1)


## Mixed state
n1<-ggplot(WEK_mixed) +
  geom_boxplot(data=WEK_mixed, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEK_mixed, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=4, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Barren state
n2<-ggplot(WEK_urchin) +
  geom_boxplot(data=WEK_urchin, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEK_urchin, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=4, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n2)$data[[1]]
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Algae state
n3<-ggplot(WEK_algae) +
  geom_boxplot(data=WEK_algae, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEK_algae, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=4, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n3)$data[[1]]
n3 <- n3 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n3 <- n3 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## locations for label plotting 
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotations
wek <- grobTree(text_grob("West End Kelp", x=.3, y=.95, hjust=0, size = 13, color = "black"))
A <- grobTree(text_grob("A", x=.25, y=.045, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.475, y=.045, hjust=0, size = 13, color = "black"))
B <- grobTree(text_grob("B", x=.725, y=.045, hjust=0, size = 13, color = "black"))


## print final plot
WEK_out <- WEK_kd_print + 
  annotation_custom(ggplotGrob(n2), xmin = (high3-width), xmax = (high3+width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n1), xmin = mid3-width, xmax = mid3+width, ymin = 0.0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n3), xmin = low3-width, xmax = low3+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  annotation_custom(wek) +
  annotation_custom(M) +
  annotation_custom(B) +
  annotation_custom(A)

print(WEK_out)
## END WestEnd Kelp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot WestEnd Urchin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inspect kernal density
WEU_kd <- ggplot(WEU, aes(NMDS1)) + geom_density(adjust = .8, color="#7b6800", size=.75,fill="lightgray",alpha=.7) +
  scale_x_continuous(limits = c(-1.5, 1.2)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme_bw() + ylab("Velocity") 
print(WEU_kd)


## empty plot for final product 
WEU_print <- ggplot(WEU, aes(NMDS1)) + 
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 


## extract data from kernal density 
p<-ggplot_build(WEU_kd)
head(p$data[[1]], 63)
kd<-as.data.frame(p$data[[1]])


## x and y from NMDS kernal density
S1 <- kd$x
S2 <- kd$y


## max y and x value from kernal density 
max_y1 <- which.max(S2)
max_x1 <- S1[max_y1]


## 2nd y and x max value 
below <- max(S2[S1 > .4])
belowY<-which(S2==below)
max_x2 <- S1[belowY]


## calculate min value between states
trough <- min(S2[S1 > max_x1 & S1 < max_x2])
ymin<-which(S2==trough)
min_x1<-S1[ymin]


## id third max ot
below2 <- max(S2[S1 < -.74])
belowY2<-which(S2==below2)
max_x3 <- S1[belowY2]


## calculate 2nd min 
trough <- min(S2[S1 > max_x3 & S1 < max_x1])
ymin2<-which(S2==trough)
min_x2<-S1[ymin2]


## check plot 
WEU_kd +
  geom_vline(xintercept=max_x1) +
  geom_vline(xintercept=max_x2) +
  geom_vline(xintercept=min_x1) + 
  geom_vline(xintercept=max_x3) +
  geom_vline(xintercept=min_x2) 


## new data frames for individual states
WEU_urchin <- filter(WEU, NMDS1>min_x1)
WEU_mixed <- filter(WEU, NMDS1<min_x1&NMDS1>min_x2)
WEU_algae <- filter(WEU, NMDS1<min_x2)


## Mixed state
n1<-ggplot(WEU_mixed) +
  geom_boxplot(data=WEU_mixed, aes(y=totalUrchin), fill="#F7F7F7",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEU_mixed, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Barren state
n2<-ggplot(WEU_urchin) +
  geom_boxplot(data=WEU_urchin, aes(y=totalUrchin), fill="#F7F7F7",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEU_urchin, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2",color="black",size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n2)$data[[1]]
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Algae state
n3<-ggplot(WEU_algae) +
  geom_boxplot(data=WEU_algae, aes(y=totalUrchin), fill="#F7F7F7",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WEU_algae, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black",size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n3)$data[[1]]
n3 <- n3 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n3 <- n3 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## point to plot labels 
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotations
weu <- grobTree(text_grob("West End Urchin", x=.3, y=.95, hjust=0, size = 13, color = "black"))
A <- grobTree(text_grob("A", x=.25, y=.045, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.475, y=.045, hjust=0, size = 13, color = "black"))
B <- grobTree(text_grob("B", x=.725, y=.045, hjust=0, size = 13, color = "black"))


## final plot 
WEU_out <- WEU_print + 
  annotation_custom(ggplotGrob(n2), xmin = (high3-width), xmax = (high3+width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n1), xmin = mid3-width, xmax = mid3+width, ymin = 0.0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n3), xmin = low3-width, xmax = low3+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  annotation_custom(weu) +
  annotation_custom(M) +
  annotation_custom(B) +
  annotation_custom(A)

print(WEU_out)
## END WestEnd Urchin plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Daytona plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## check kernal density
Day_kd <- ggplot(Day, aes(NMDS1)) +
  geom_density()
print(Day_kd)


## setup final plot
Day_kd_print <- ggplot(Day, aes(NMDS1)) +
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 


## extract data from kernal density 
p<-ggplot_build(Day_kd)
head(p$data[[1]], 63)
kd<-as.data.frame(p$data[[1]])


### x and y of NMDS kernal density
S1 <- kd$x; 
S2 <- kd$y


## x and y max value from kernal density 
max_y1 <- which.max(S2); 
max_x1 <- S1[max_y1]


## id second max 
below <- max(S2[S1>.5]);
belowY<-which(S2==below); 
max_x2 <- S1[belowY]


## calculate min 
trough <- min(S2[S1 > max_x1 & S1 < max_x2]); 
ymin<-which(S2==trough); 
min_x1<-S1[ymin]


## check plot 
ggplot(Day, aes(NMDS1)) + geom_density() + 
  geom_vline(xintercept=max_x1) + 
  geom_vline(xintercept=max_x2) +   
  geom_vline(xintercept=min_x1)


## new data frames
Day_urchin <- filter(Day, NMDS1>min_x1); 
Day_mixed <- filter(Day, NMDS1<min_x1)


## Mixed state
n1<-ggplot(Day_mixed) +
  geom_boxplot(data=Day_mixed, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,2600) + 
  geom_beeswarm(data=Day_mixed, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Barren state
n2<-ggplot(Day_urchin) +
  geom_boxplot(data=Day_urchin, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,2600) + 
  geom_beeswarm(data=Day_urchin, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n2)$data[[1]]
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## custom label positioning 
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotations 
dy <- grobTree(text_grob("Daytona", x=.4, y=.95, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.33333, y=.045, hjust=0, size = 13, color = "black"))
B <- grobTree(text_grob("B", x=.63333, y=.045, hjust=0, size = 13, color = "black"))


## final plot 
Day_out<-Day_kd_print + 
  annotation_custom(ggplotGrob(n2), xmin = (high2-width), xmax = (high2+width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n1), xmin = low2-width, xmax = low2+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  annotation_custom(dy) +
  annotation_custom(M) +
  annotation_custom(B) 

print(Day_out)
## END Daytona plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## East Dutch plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## inspect kernal density
ED_kd <- ggplot(ED, aes(NMDS1)) +
  geom_density()
print(ED_kd)


## setup final plot 
ED_kd_print <- ggplot(ED, aes(NMDS1)) + 
  scale_x_continuous(limits = c(0,1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 


## extract data from kernal density 
p<-ggplot_build(ED_kd)
head(p$data[[1]], 63)
kd<-as.data.frame(p$data[[1]])


### x and y of NMDS kernal density
S1 <- kd$x; 
S2 <- kd$y


## maX value of kernal y; corresponds to what value X
max_y1 <- which.max(S2); 
max_x1 <- S1[max_y1]


## id second max 
below <- max(S2[S1 < -0.5]);
belowY<-which(S2==below); 
max_x2 <- S1[belowY]


## calculate min 
trough <- min(S2[S1 < max_x1 & S1 > max_x2]); 
ymin<-which(S2==trough); 
min_x1<-S1[ymin]


## check plot 
ggplot(ED, aes(NMDS1)) + geom_density() + 
  geom_vline(xintercept=max_x1) + 
  geom_vline(xintercept=max_x2) +   
  geom_vline(xintercept=min_x1)


## data frames
ED_mixed <- filter(ED, NMDS1>min_x1); 
ED_algae <- filter(ED, NMDS1<min_x1)


## Mixed state
n1<-ggplot(ED) +
  geom_boxplot(data=ED_mixed, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=ED_mixed, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## Algae state
n2<-ggplot(ED) +
  geom_boxplot(data=ED_algae, aes(y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=ED_algae, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n2)$data[[1]]
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n2 <- n2 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## setup annotation locations 
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotations 
ed <- grobTree(text_grob("East Dutch", x=.35, y=.95, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.33333, y=.045, hjust=0, size = 13, color = "black"))
B <- grobTree(text_grob("B", x=.63333, y=.045, hjust=0, size = 13, color = "black"))


## final plot
ED_out <- ED_kd_print + 
  annotation_custom(ggplotGrob(n1), xmin = high2-width, xmax = high2+width, ymin = 0.0, ymax = 2.5) +
  annotation_custom(ggplotGrob(n2), xmin = low2-width, xmax = low2+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  annotation_custom(ed) +
  annotation_custom(M) +
  annotation_custom(B) 

print(ED_out)
## END East Dutch plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot West Dutch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
WD_kd <- ggplot(WD, aes(NMDS1)) +
  geom_density()
print(WD_kd)


WD_kd_print <- ggplot(WD, aes(NMDS1)) + 
  scale_x_continuous(limits = c(0, 1)) + scale_y_continuous(limits = c(0,2.5), expand = c(0,0)) + theme_bw() + ylab("Velocity") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        plot.title = element_blank()) + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) + theme(legend.position="none") 


## single, Mixed state
n1<-ggplot(WD) +
  geom_boxplot(data=WD, aes(x=0,y=totalUrchin), fill="white",alpha=0, color="white", outlier.shape=NA) + ylim(0,1600) + 
  geom_beeswarm(data=WD, aes(x=0, y=totalUrchin), shape=21, fill ="#C2C2C2", color="black", size=.8,stroke=.35,
                groupOnX = NULL, cex=5, priority = c("random")) +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA)) +
  theme(plot.title = element_blank(),axis.title.y = element_blank(),axis.title.x = element_blank(),axis.text.y.left = element_blank(),axis.text.x = element_blank(),axis.ticks = element_blank(),legend.position = "none") 
fig_dat <- ggplot_build(n1)$data[[1]]
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1.1, xend=0+1.1, y=middle, yend=middle), colour="black", size=1.5)
n1 <- n1 + geom_segment(data=fig_dat, aes(x=0-1, xend=0+1, y=middle, yend=middle), colour="red", size=1)


## custom annotation locations 
plot_xmin <- 0; plot_xmax <- 1; low2 <- 0.33333; high2 <- 0.66666; low3<-0.25; mid3<-0.5; high3<-0.75;
width<-(plot_xmax-plot_xmin)/5.75


## custom annotations 
wd <- grobTree(text_grob("West Dutch", x=.35, y=.95, hjust=0, size = 13, color = "black"))
M <- grobTree(text_grob("M", x=.48, y=.045, hjust=0, size = 13, color = "black"))


## final plot
WD_out <- WD_kd_print + 
  annotation_custom(ggplotGrob(n1), xmin = mid3-width, xmax = mid3+width, ymin = 0.0, ymax = 2.5) +
  theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) +
  annotation_custom(wd) +
  annotation_custom(M) 

print(WD_out)
## END West Dutch plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Aggregate all plots for final ms figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=3,w=15, record=TRUE)


## custom axis for later call 
new<-grid.yaxis(at=c(0,400,800,1200,1600), 
                vp=vpStack(viewport(width=unit(2,"lines")),
                           viewport(x=1, yscale = c(-50,1050), just="left")))


## aggregate all plot 
p2 <- ggarrange(tag_facet(NF_out +
                            facet_wrap(~"NMDS1"), 
                          tag_pool = "a"),
                tag_facet(WEK_out +
                            facet_wrap(~"NMDS1"),
                          tag_pool = "b"),
                tag_facet(WEU_out +
                            facet_wrap(~"NMDS1"),
                          tag_pool = "c"),
                tag_facet(Day_out +
                            facet_wrap(~"NMDS1"),
                          tag_pool = "d"),
                tag_facet(ED_out +
                            facet_wrap(~"NMDS1"),
                          tag_pool = "e" ),
                tag_facet(WD_out +
                            facet_wrap(~"NMDS1"),
                          tag_pool = "f" ),
                nrow=1)


## add axis legend 
name<-grid.text(label="total Urchin abundance", x=.01,y=.5,rot=90)

## add previously created axis
sample_vp <- viewport(x = .545, y = 0.3475, width = 1, height = .575, just = c("right", "center"))
pushViewport(sample_vp)
grid.draw(new)
grid.draw(name)
popViewport()   
## END OF SCRIPT ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


