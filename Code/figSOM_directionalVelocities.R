## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Figure for Substrate_Complexity SOM figure: Directional velocities of community shift
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## Initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
#setwd("D:/SNI/Data/SEP2018")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16)) 

vel.theme = theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  panel.background = element_rect(fill = "transparent",colour = NA), 
                  plot.background = element_rect(fill = "transparent",colour = NA), 
                  plot.title = element_blank(), 
                  axis.title.y = element_blank(), 
                  axis.title.x = element_blank(), 
                  axis.text.y.left = element_blank(), 
                  axis.text.x = element_blank(), 
                  axis.ticks = element_blank()) 

no.axes = theme(axis.title = element_blank(),
                axis.text = element_blank(),
                axis.ticks = element_blank())

show.left = theme(axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.ticks.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks.y = element_blank())

show.below = theme(axis.title.y = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks.y = element_blank(),
                   axis.text.x = element_blank(),
                   axis.ticks.x = element_blank())
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## data wrangling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## directional-specific velocities
dat$left <- ifelse(dat$X_Vel<0, dat$X_Vel, NA)
dat$right <- ifelse(dat$X_Vel>0, dat$X_Vel, NA)


## midpts to match directional velocities 
dat$leftMid <- ifelse(dat$X_Vel<0, dat$MidPt, NA)
dat$rightMid <- ifelse(dat$X_Vel>0, dat$MidPt, NA)


## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))


## specific transects
NF_10R <- filter(NF, TRANSECT %in% c("10R")) 
NF_22R <- filter(NF, TRANSECT %in% c("22R")) 
NF_32L <- filter(NF, TRANSECT %in% c("32L")) 
NF_39R <- filter(NF, TRANSECT %in% c("39R")) 
NF_45R <- filter(NF, TRANSECT %in% c("45R")) 

WEK_10R <- filter(WEK, TRANSECT %in% c("10R")) 
WEK_22R <- filter(WEK, TRANSECT %in% c("22R")) 
WEK_32L <- filter(WEK, TRANSECT %in% c("32L")) 
WEK_39R <- filter(WEK, TRANSECT %in% c("39R")) 
WEK_45L <- filter(WEK, TRANSECT %in% c("45L")) 

WEU_10L <- filter(WEU, TRANSECT %in% c("10L")) 
WEU_22L <- filter(WEU, TRANSECT %in% c("22L")) 
WEU_32R <- filter(WEU, TRANSECT %in% c("32R")) 
WEU_39L <- filter(WEU, TRANSECT %in% c("39L")) 
WEU_45L <- filter(WEU, TRANSECT %in% c("45L")) 

WD_10R <- filter(WD, TRANSECT %in% c("10R"))
WD_22L <- filter(WD, TRANSECT %in% c("22L"))
WD_32L <- filter(WD, TRANSECT %in% c("32L"))
WD_39L <- filter(WD, TRANSECT %in% c("39L"))
WD_45L <- filter(WD, TRANSECT %in% c("45L"))

ED_10R <- filter(ED, TRANSECT %in% c("10R"))
ED_22R <- filter(ED, TRANSECT %in% c("22R"))
ED_32L <- filter(ED, TRANSECT %in% c("32L"))
ED_39R <- filter(ED, TRANSECT %in% c("39R"))
ED_45R <- filter(ED, TRANSECT %in% c("45R"))

Day_10R <- filter(Day, TRANSECT %in% c("10R"))
Day_22L <- filter(Day, TRANSECT %in% c("22L"))
Day_22R <- filter(Day, TRANSECT %in% c("22R"))
Day_32L <- filter(Day, TRANSECT %in% c("32L"))
Day_39L <- filter(Day, TRANSECT %in% c("39L"))
## END subsetting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## universal graphing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## all subfigs
xmin <- -1.6
xmax <- 1.6

## kernal density options 
fillcolor <- "#CDCDCD"
lineWeight_kd <- .75
ymin_kd <- 0
ymax_kd <- 3
ymax_kd_day <- 3
kernalAdjust <- 1

## loess options
loessSpan <- 0.75
loessLeft <- "#008000"
loessRight <- "#660198"
lineWeight_loess <- .75
ymin_line <- 0
ymax_line <- 0.01
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## transect specific colors ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NF_10R_col <- "#390b0b"
NF_22R_col <- "#721716"
NF_32L_col <- "#ab2221"
NF_39R_col <- "#cb5151"
NF_45R_col <- "#df9392"

WEK_10R_col <- "#301800"  
WEK_22R_col <- "#773b00"
WEK_32L_col <- "#be5e00"
WEK_39R_col <- "#f0841a"
WEK_45L_col <- "#f5ad66"

WEU_10L_col <- "#141100"  
WEU_22L_col <- "#3d3400"
WEU_32R_col <- "#7b6800"
WEU_39L_col <- "#b99c00"
WEU_45L_col <- "#e6d680"

Day_10R_col <- "#131a0d"  
Day_22L_col <- "#304221"
Day_22R_col <- "#4d6934"
Day_32L_col <- "#708f54"
Day_39L_col <- "#a0b58d"

ED_10R_col <- "#003446" 
ED_22R_col <- "#00536f"
ED_32L_col <- "#3386a2"
ED_39R_col <- "#66a4b9"
ED_45R_col <- "#99c3d1"

WD_10R_col <- "#2b112b"  
WD_22L_col <- "#552255"
WD_32L_col <- "#803280"
WD_39L_col <- "#a560a5"
WD_45L_col <- "#c79cc7"
## END graphing options ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







graphics.off()
windows(h=4,w=4, record = TRUE)



## NavFac ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
NF_10R_left <- filter(NF_10R, left < 0) 
NF_10R_right <- filter(NF_10R, right > 0)
## velocities
smooth_NF_10R_left <- loess(NF_10R_left$Vel ~ NF_10R_left$leftMid, data = NF_10R, span=loessSpan)
loess_NF_10R_left <- predict(smooth_NF_10R_left)
smooth_NF_10R_right <- loess(NF_10R_right$Vel ~ NF_10R_right$rightMid, data = NF_10R, span=loessSpan)
loess_NF_10R_right <- predict(smooth_NF_10R_right)
## kernal density 
NF_kd_10R <- ggplot(NF_10R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=NF_10R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
NF_10R_kelp <- ggplot(NF_10R_left, aes(x=MidPt, y=loess_NF_10R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
NF_10R_urchin <- ggplot(NF_10R_right, aes(x=MidPt, y=loess_NF_10R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
NF_10R_both <- NF_kd_10R + annotation_custom(ggplotGrob(NF_10R_kelp)) + annotation_custom(ggplotGrob(NF_10R_urchin))
print(NF_10R_both)
  
  
## subset based on direction
NF_22R_left <- filter(NF_22R, left < 0) 
NF_22R_right <- filter(NF_22R, right > 0)
## velocities
smooth_NF_22R_left <- loess(NF_22R_left$Vel ~ NF_22R_left$leftMid, data = NF_22R, span=loessSpan)
loess_NF_22R_left <- predict(smooth_NF_22R_left)
smooth_NF_22R_right <- loess(NF_22R_right$Vel ~ NF_22R_right$rightMid, data = NF_22R, span=loessSpan)
loess_NF_22R_right <- predict(smooth_NF_22R_right)
## kernal density 
NF_kd_22R <- ggplot(NF_22R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=NF_22R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
NF_22R_kelp <- ggplot(NF_22R_left, aes(x=MidPt, y=loess_NF_22R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
NF_22R_urchin <- ggplot(NF_22R_right, aes(x=MidPt, y=loess_NF_22R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
NF_22R_both <- NF_kd_22R + annotation_custom(ggplotGrob(NF_22R_kelp)) + annotation_custom(ggplotGrob(NF_22R_urchin))
print(NF_22R_both)


## subset based on direction
NF_32L_left <- filter(NF_32L, left < 0) 
NF_32L_right <- filter(NF_32L, right > 0)
## velocities
smooth_NF_32L_left <- loess(NF_32L_left$Vel ~ NF_32L_left$leftMid, data = NF_32L, span=loessSpan)
loess_NF_32L_left <- predict(smooth_NF_32L_left)
smooth_NF_32L_right <- loess(NF_32L_right$Vel ~ NF_32L_right$rightMid, data = NF_32L, span=loessSpan)
loess_NF_32L_right <- predict(smooth_NF_32L_right)
## kernal density 
NF_kd_32L <- ggplot(NF_32L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=NF_32L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
NF_32L_kelp <- ggplot(NF_32L_left, aes(x=MidPt, y=loess_NF_32L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
NF_32L_urchin <- ggplot(NF_32L_right, aes(x=MidPt, y=loess_NF_32L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
NF_32L_both <- NF_kd_32L + annotation_custom(ggplotGrob(NF_32L_kelp)) + annotation_custom(ggplotGrob(NF_32L_urchin))
print(NF_32L_both)


## subset based on direction
NF_39R_left <- filter(NF_39R, left < 0) 
NF_39R_right <- filter(NF_39R, right > 0)
## velocities
smooth_NF_39R_left <- loess(NF_39R_left$Vel ~ NF_39R_left$leftMid, data = NF_39R, span=loessSpan)
loess_NF_39R_left <- predict(smooth_NF_39R_left)
smooth_NF_39R_right <- loess(NF_39R_right$Vel ~ NF_39R_right$rightMid, data = NF_39R, span=loessSpan)
loess_NF_39R_right <- predict(smooth_NF_39R_right)
## kernal density 
NF_kd_39R <- ggplot(NF_39R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=NF_39R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
NF_39R_kelp <- ggplot(NF_39R_left, aes(x=MidPt, y=loess_NF_39R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
NF_39R_urchin <- ggplot(NF_39R_right, aes(x=MidPt, y=loess_NF_39R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
NF_39R_both <- NF_kd_39R + annotation_custom(ggplotGrob(NF_39R_kelp)) + annotation_custom(ggplotGrob(NF_39R_urchin))
print(NF_39R_both)


## subset based on direction
NF_45R_left <- filter(NF_45R, left < 0) 
NF_45R_right <- filter(NF_45R, right > 0)
## velocities
smooth_NF_45R_left <- loess(NF_45R_left$Vel ~ NF_45R_left$leftMid, data = NF_45R, span=loessSpan)
loess_NF_45R_left <- predict(smooth_NF_45R_left)
smooth_NF_45R_right <- loess(NF_45R_right$Vel ~ NF_45R_right$rightMid, data = NF_45R, span=loessSpan)
loess_NF_45R_right <- predict(smooth_NF_45R_right)
## kernal density 
NF1 <- grobTree(text_grob("NavFac", x=.35, y=.75, hjust=0, size = 12, color = "black"))
NF_kd_45R <- ggplot(NF_45R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=NF_45R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes +
  annotation_custom(NF1)
## plot velocity lines
NF_45R_kelp <- ggplot(NF_45R_left, aes(x=MidPt, y=loess_NF_45R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
NF_45R_urchin <- ggplot(NF_45R_right, aes(x=MidPt, y=loess_NF_45R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
NF_45R_both <- NF_kd_45R + annotation_custom(ggplotGrob(NF_45R_kelp)) + annotation_custom(ggplotGrob(NF_45R_urchin))
print(NF_45R_both)
## END NavFac ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## WestEnd Kelp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
WEK_10R_left <- filter(WEK_10R, left < 0) 
WEK_10R_right <- filter(WEK_10R, right > 0)
## velocities
smooth_WEK_10R_left <- loess(WEK_10R_left$Vel ~ WEK_10R_left$leftMid, data = WEK_10R, span=loessSpan)
loess_WEK_10R_left <- predict(smooth_WEK_10R_left)
smooth_WEK_10R_right <- loess(WEK_10R_right$Vel ~ WEK_10R_right$rightMid, data = WEK_10R, span=loessSpan)
loess_WEK_10R_right <- predict(smooth_WEK_10R_right)
## kernal density 
WEK_kd_10R <- ggplot(WEK_10R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEK_10R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEK_10R_kelp <- ggplot(WEK_10R_left, aes(x=MidPt, y=loess_WEK_10R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEK_10R_urchin <- ggplot(WEK_10R_right, aes(x=MidPt, y=loess_WEK_10R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEK_10R_both <- WEK_kd_10R + annotation_custom(ggplotGrob(WEK_10R_kelp)) + annotation_custom(ggplotGrob(WEK_10R_urchin))
print(WEK_10R_both)


## subset based on direction
WEK_22R_left <- filter(WEK_22R, left < 0) 
WEK_22R_right <- filter(WEK_22R, right > 0)
## velocities
smooth_WEK_22R_left <- loess(WEK_22R_left$Vel ~ WEK_22R_left$leftMid, data = WEK_22R, span=loessSpan)
loess_WEK_22R_left <- predict(smooth_WEK_22R_left)
smooth_WEK_22R_right <- loess(WEK_22R_right$Vel ~ WEK_22R_right$rightMid, data = WEK_22R, span=loessSpan)
loess_WEK_22R_right <- predict(smooth_WEK_22R_right)
## kernal density 
WEK_kd_22R <- ggplot(WEK_22R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEK_22R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEK_22R_kelp <- ggplot(WEK_22R_left, aes(x=MidPt, y=loess_WEK_22R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEK_22R_urchin <- ggplot(WEK_22R_right, aes(x=MidPt, y=loess_WEK_22R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEK_22R_both <- WEK_kd_22R + annotation_custom(ggplotGrob(WEK_22R_kelp)) + annotation_custom(ggplotGrob(WEK_22R_urchin))
print(WEK_22R_both)


## subset based on direction
WEK_32L_left <- filter(WEK_32L, left < 0) 
WEK_32L_right <- filter(WEK_32L, right > 0)
## velocities
smooth_WEK_32L_left <- loess(WEK_32L_left$Vel ~ WEK_32L_left$leftMid, data = WEK_32L, span=loessSpan)
loess_WEK_32L_left <- predict(smooth_WEK_32L_left)
smooth_WEK_32L_right <- loess(WEK_32L_right$Vel ~ WEK_32L_right$rightMid, data = WEK_32L, span=loessSpan)
loess_WEK_32L_right <- predict(smooth_WEK_32L_right)
## kernal density 
WEK_kd_32L <- ggplot(WEK_32L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEK_32L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEK_32L_kelp <- ggplot(WEK_32L_left, aes(x=MidPt, y=loess_WEK_32L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEK_32L_urchin <- ggplot(WEK_32L_right, aes(x=MidPt, y=loess_WEK_32L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEK_32L_both <- WEK_kd_32L + annotation_custom(ggplotGrob(WEK_32L_kelp)) + annotation_custom(ggplotGrob(WEK_32L_urchin))
print(WEK_32L_both)


## subset based on direction
WEK_39R_left <- filter(WEK_39R, left < 0) 
WEK_39R_right <- filter(WEK_39R, right > 0)
## velocities
smooth_WEK_39R_left <- loess(WEK_39R_left$Vel ~ WEK_39R_left$leftMid, data = WEK_39R, span=loessSpan)
loess_WEK_39R_left <- predict(smooth_WEK_39R_left)
smooth_WEK_39R_right <- loess(WEK_39R_right$Vel ~ WEK_39R_right$rightMid, data = WEK_39R, span=loessSpan)
loess_WEK_39R_right <- predict(smooth_WEK_39R_right)
## kernal density 
WEK_kd_39R <- ggplot(WEK_39R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEK_39R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEK_39R_kelp <- ggplot(WEK_39R_left, aes(x=MidPt, y=loess_WEK_39R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEK_39R_urchin <- ggplot(WEK_39R_right, aes(x=MidPt, y=loess_WEK_39R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEK_39R_both <- WEK_kd_39R + annotation_custom(ggplotGrob(WEK_39R_kelp)) + annotation_custom(ggplotGrob(WEK_39R_urchin))
print(WEK_39R_both)


## subset based on direction
WEK_45L_left <- filter(WEK_45L, left < 0) 
WEK_45L_right <- filter(WEK_45L, right > 0)
## velocities
smooth_WEK_45L_left <- loess(WEK_45L_left$Vel ~ WEK_45L_left$leftMid, data = WEK_45L, span=loessSpan)
loess_WEK_45L_left <- predict(smooth_WEK_45L_left)
smooth_WEK_45L_right <- loess(WEK_45L_right$Vel ~ WEK_45L_right$rightMid, data = WEK_45L, span=loessSpan)
loess_WEK_45L_right <- predict(smooth_WEK_45L_right)
## kernal density 
WEK1 <- grobTree(text_grob("WestEnd Kelp", x=.15, y=.75, hjust=0, size = 12, color = "black"))
WEK_kd_45L <- ggplot(WEK_45L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEK_45L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes +
  annotation_custom(WEK1)

## plot velocity lines
WEK_45L_kelp <- ggplot(WEK_45L_left, aes(x=MidPt, y=loess_WEK_45L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEK_45L_urchin <- ggplot(WEK_45L_right, aes(x=MidPt, y=loess_WEK_45L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEK_45L_both <- WEK_kd_45L + annotation_custom(ggplotGrob(WEK_45L_kelp)) + annotation_custom(ggplotGrob(WEK_45L_urchin))
print(WEK_45L_both)
## END WestEnd Kelp ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## WestEnd Urchin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
WEU_10L_left <- filter(WEU_10L, left < 0) 
WEU_10L_right <- filter(WEU_10L, right > 0)
## velocities
smooth_WEU_10L_left <- loess(WEU_10L_left$Vel ~ WEU_10L_left$leftMid, data = WEU_10L, span=loessSpan)
loess_WEU_10L_left <- predict(smooth_WEU_10L_left)
smooth_WEU_10L_right <- loess(WEU_10L_right$Vel ~ WEU_10L_right$rightMid, data = WEU_10L, span=loessSpan)
loess_WEU_10L_right <- predict(smooth_WEU_10L_right)
## kernal density 
WEU_kd_10L <- ggplot(WEU_10L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEU_10L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + 
  show.left + ylab("Negative potential and velocity")
## plot velocity lines
WEU_10L_kelp <- ggplot(WEU_10L_left, aes(x=MidPt, y=loess_WEU_10L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEU_10L_urchin <- ggplot(WEU_10L_right, aes(x=MidPt, y=loess_WEU_10L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEU_10L_both <- WEU_kd_10L + annotation_custom(ggplotGrob(WEU_10L_kelp)) + annotation_custom(ggplotGrob(WEU_10L_urchin))
print(WEU_10L_both)


## subset based on direction
WEU_22L_left <- filter(WEU_22L, left < 0) 
WEU_22L_right <- filter(WEU_22L, right > 0)
## velocities
smooth_WEU_22L_left <- loess(WEU_22L_left$Vel ~ WEU_22L_left$leftMid, data = WEU_22L, span=loessSpan)
loess_WEU_22L_left <- predict(smooth_WEU_22L_left)
smooth_WEU_22L_right <- loess(WEU_22L_right$Vel ~ WEU_22L_right$rightMid, data = WEU_22L, span=loessSpan)
loess_WEU_22L_right <- predict(smooth_WEU_22L_right)
## kernal density 
WEU_kd_22L <- ggplot(WEU_22L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEU_22L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEU_22L_kelp <- ggplot(WEU_22L_left, aes(x=MidPt, y=loess_WEU_22L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEU_22L_urchin <- ggplot(WEU_22L_right, aes(x=MidPt, y=loess_WEU_22L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEU_22L_both <- WEU_kd_22L + annotation_custom(ggplotGrob(WEU_22L_kelp)) + annotation_custom(ggplotGrob(WEU_22L_urchin))
print(WEU_22L_both)


## subset based on direction
WEU_32R_left <- filter(WEU_32R, left < 0) 
WEU_32R_right <- filter(WEU_32R, right > 0)
## velocities
smooth_WEU_32R_left <- loess(WEU_32R_left$Vel ~ WEU_32R_left$leftMid, data = WEU_32R, span=loessSpan)
loess_WEU_32R_left <- predict(smooth_WEU_32R_left)
smooth_WEU_32R_right <- loess(WEU_32R_right$Vel ~ WEU_32R_right$rightMid, data = WEU_32R, span=loessSpan)
loess_WEU_32R_right <- predict(smooth_WEU_32R_right)
## kernal density 
WEU_kd_32R <- ggplot(WEU_32R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEU_32R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEU_32R_kelp <- ggplot(WEU_32R_left, aes(x=MidPt, y=loess_WEU_32R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEU_32R_urchin <- ggplot(WEU_32R_right, aes(x=MidPt, y=loess_WEU_32R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEU_32R_both <- WEU_kd_32R + annotation_custom(ggplotGrob(WEU_32R_kelp)) + annotation_custom(ggplotGrob(WEU_32R_urchin))
print(WEU_32R_both)


## subset based on direction
WEU_39L_left <- filter(WEU_39L, left < 0) 
WEU_39L_right <- filter(WEU_39L, right > 0)
## velocities
smooth_WEU_39L_left <- loess(WEU_39L_left$Vel ~ WEU_39L_left$leftMid, data = WEU_39L, span=loessSpan)
loess_WEU_39L_left <- predict(smooth_WEU_39L_left)
smooth_WEU_39L_right <- loess(WEU_39L_right$Vel ~ WEU_39L_right$rightMid, data = WEU_39L, span=loessSpan)
loess_WEU_39L_right <- predict(smooth_WEU_39L_right)
## kernal density 
WEU_kd_39L <- ggplot(WEU_39L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEU_39L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes
## plot velocity lines
WEU_39L_kelp <- ggplot(WEU_39L_left, aes(x=MidPt, y=loess_WEU_39L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEU_39L_urchin <- ggplot(WEU_39L_right, aes(x=MidPt, y=loess_WEU_39L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEU_39L_both <- WEU_kd_39L + annotation_custom(ggplotGrob(WEU_39L_kelp)) + annotation_custom(ggplotGrob(WEU_39L_urchin))
print(WEU_39L_both)


## subset based on direction
WEU_45L_left <- filter(WEU_45L, left < 0) 
WEU_45L_right <- filter(WEU_45L, right > 0)
## velocities
smooth_WEU_45L_left <- loess(WEU_45L_left$Vel ~ WEU_45L_left$leftMid, data = WEU_45L, span=loessSpan)
loess_WEU_45L_left <- predict(smooth_WEU_45L_left)
smooth_WEU_45L_right <- loess(WEU_45L_right$Vel ~ WEU_45L_right$rightMid, data = WEU_45L, span=loessSpan)
loess_WEU_45L_right <- predict(smooth_WEU_45L_right)
## kernal density 
WEU1 <- grobTree(text_grob("WestEnd Urchin", x=.10, y=.75, hjust=0, size = 12, color = "black"))
WEU_kd_45L <- ggplot(WEU_45L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WEU_45L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd)) + no.axes +
  annotation_custom(WEU1)
## plot velocity lines
WEU_45L_kelp <- ggplot(WEU_45L_left, aes(x=MidPt, y=loess_WEU_45L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WEU_45L_urchin <- ggplot(WEU_45L_right, aes(x=MidPt, y=loess_WEU_45L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WEU_45L_both <- WEU_kd_45L + annotation_custom(ggplotGrob(WEU_45L_kelp)) + annotation_custom(ggplotGrob(WEU_45L_urchin))
print(WEU_45L_both)
## END WestEnd Urchin ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Daytona ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
Day_10R_left <- filter(Day_10R, left < 0) 
Day_10R_right <- filter(Day_10R, right > 0)
## velocities
smooth_Day_10R_left <- loess(Day_10R_left$Vel ~ Day_10R_left$leftMid, data = Day_10R, span=loessSpan)
loess_Day_10R_left <- predict(smooth_Day_10R_left)
smooth_Day_10R_right <- loess(Day_10R_right$Vel ~ Day_10R_right$rightMid, data = Day_10R, span=loessSpan)
loess_Day_10R_right <- predict(smooth_Day_10R_right)
## kernal density 
Day_kd_10R <- ggplot(Day_10R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=Day_10R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
Day_10R_kelp <- ggplot(Day_10R_left, aes(x=MidPt, y=loess_Day_10R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
Day_10R_urchin <- ggplot(Day_10R_right, aes(x=MidPt, y=loess_Day_10R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
Day_10R_both <- Day_kd_10R + annotation_custom(ggplotGrob(Day_10R_kelp)) + annotation_custom(ggplotGrob(Day_10R_urchin))
print(Day_10R_both)


## subset based on direction
Day_22L_left <- filter(Day_22L, left < 0) 
Day_22L_right <- filter(Day_22L, right > 0)
## velocities
smooth_Day_22L_left <- loess(Day_22L_left$Vel ~ Day_22L_left$leftMid, data = Day_22L, span=loessSpan)
loess_Day_22L_left <- predict(smooth_Day_22L_left)
smooth_Day_22L_right <- loess(Day_22L_right$Vel ~ Day_22L_right$rightMid, data = Day_22L, span=loessSpan)
loess_Day_22L_right <- predict(smooth_Day_22L_right)
## kernal density 
Day_kd_22L <- ggplot(Day_22L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=Day_22L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
Day_22L_kelp <- ggplot(Day_22L_left, aes(x=MidPt, y=loess_Day_22L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
Day_22L_urchin <- ggplot(Day_22L_right, aes(x=MidPt, y=loess_Day_22L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
Day_22L_both <- Day_kd_22L + annotation_custom(ggplotGrob(Day_22L_kelp)) + annotation_custom(ggplotGrob(Day_22L_urchin))
print(Day_22L_both)


## subset based on direction
Day_22R_left <- filter(Day_22R, left < 0) 
Day_22R_right <- filter(Day_22R, right > 0)
## velocities
smooth_Day_22R_left <- loess(Day_22R_left$Vel ~ Day_22R_left$leftMid, data = Day_22R, span=loessSpan)
loess_Day_22R_left <- predict(smooth_Day_22R_left)
smooth_Day_22R_right <- loess(Day_22R_right$Vel ~ Day_22R_right$rightMid, data = Day_22R, span=loessSpan)
loess_Day_22R_right <- predict(smooth_Day_22R_right)
## kernal density 
Day_kd_22R <- ggplot(Day_22R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=Day_22R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
Day_22R_kelp <- ggplot(Day_22R_left, aes(x=MidPt, y=loess_Day_22R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
Day_22R_urchin <- ggplot(Day_22R_right, aes(x=MidPt, y=loess_Day_22R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
Day_22R_both <- Day_kd_22R + annotation_custom(ggplotGrob(Day_22R_kelp)) + annotation_custom(ggplotGrob(Day_22R_urchin))
print(Day_22R_both)


## subset based on direction
Day_32L_left <- filter(Day_32L, left < 0) 
Day_32L_right <- filter(Day_32L, right > 0)
## velocities
smooth_Day_32L_left <- loess(Day_32L_left$Vel ~ Day_32L_left$leftMid, data = Day_32L, span=loessSpan)
loess_Day_32L_left <- predict(smooth_Day_32L_left)
smooth_Day_32L_right <- loess(Day_32L_right$Vel ~ Day_32L_right$rightMid, data = Day_32L, span=loessSpan)
loess_Day_32L_right <- predict(smooth_Day_32L_right)
## kernal density 
Day_kd_32L <- ggplot(Day_32L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=Day_32L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
Day_32L_kelp <- ggplot(Day_32L_left, aes(x=MidPt, y=loess_Day_32L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
Day_32L_urchin <- ggplot(Day_32L_right, aes(x=MidPt, y=loess_Day_32L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
Day_32L_both <- Day_kd_32L + annotation_custom(ggplotGrob(Day_32L_kelp)) + annotation_custom(ggplotGrob(Day_32L_urchin))
print(Day_32L_both)


## subset based on direction
Day_39L_left <- filter(Day_39L, left < 0) 
Day_39L_right <- filter(Day_39L, right > 0)
## velocities
smooth_Day_39L_left <- loess(Day_39L_left$Vel ~ Day_39L_left$leftMid, data = Day_39L, span=loessSpan)
loess_Day_39L_left <- predict(smooth_Day_39L_left)
smooth_Day_39L_right <- loess(Day_39L_right$Vel ~ Day_39L_right$rightMid, data = Day_39L, span=loessSpan)
loess_Day_39L_right <- predict(smooth_Day_39L_right)
## kernal density 
Day1 <- grobTree(text_grob("Daytona", x=.35, y=.75, hjust=0, size = 12, color = "black"))
Day_kd_39L <- ggplot(Day_39L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=Day_39L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes +
  annotation_custom(Day1)

## plot velocity lines
Day_39L_kelp <- ggplot(Day_39L_left, aes(x=MidPt, y=loess_Day_39L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
Day_39L_urchin <- ggplot(Day_39L_right, aes(x=MidPt, y=loess_Day_39L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
Day_39L_both <- Day_kd_39L + annotation_custom(ggplotGrob(Day_39L_kelp)) + annotation_custom(ggplotGrob(Day_39L_urchin))
print(Day_39L_both)
## END Daytona ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## East Dutch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
ED_10R_left <- filter(ED_10R, left < 0) 
ED_10R_right <- filter(ED_10R, right > 0)
## velocities
smooth_ED_10R_left <- loess(ED_10R_left$Vel ~ ED_10R_left$leftMid, data = ED_10R, span=loessSpan)
loess_ED_10R_left <- predict(smooth_ED_10R_left)
smooth_ED_10R_right <- loess(ED_10R_right$Vel ~ ED_10R_right$rightMid, data = ED_10R, span=loessSpan)
loess_ED_10R_right <- predict(smooth_ED_10R_right)
## kernal density 
ED_kd_10R <- ggplot(ED_10R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=ED_10R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
ED_10R_kelp <- ggplot(ED_10R_left, aes(x=MidPt, y=loess_ED_10R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
ED_10R_urchin <- ggplot(ED_10R_right, aes(x=MidPt, y=loess_ED_10R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
ED_10R_both <- ED_kd_10R + annotation_custom(ggplotGrob(ED_10R_kelp)) + annotation_custom(ggplotGrob(ED_10R_urchin))
print(ED_10R_both)


## subset based on direction
ED_22R_left <- filter(ED_22R, left < 0) 
ED_22R_right <- filter(ED_22R, right > 0)
## velocities
smooth_ED_22R_left <- loess(ED_22R_left$Vel ~ ED_22R_left$leftMid, data = ED_22R, span=loessSpan)
loess_ED_22R_left <- predict(smooth_ED_22R_left)
smooth_ED_22R_right <- loess(ED_22R_right$Vel ~ ED_22R_right$rightMid, data = ED_22R, span=loessSpan)
loess_ED_22R_right <- predict(smooth_ED_22R_right)
## kernal density 
ED_kd_22R <- ggplot(ED_22R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=ED_22R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
ED_22R_kelp <- ggplot(ED_22R_left, aes(x=MidPt, y=loess_ED_22R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
ED_22R_urchin <- ggplot(ED_22R_right, aes(x=MidPt, y=loess_ED_22R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
ED_22R_both <- ED_kd_22R + annotation_custom(ggplotGrob(ED_22R_kelp)) + annotation_custom(ggplotGrob(ED_22R_urchin))
print(ED_22R_both)


## subset based on direction
ED_32L_left <- filter(ED_32L, left < 0) 
ED_32L_right <- filter(ED_32L, right > 0)
## velocities
smooth_ED_32L_left <- loess(ED_32L_left$Vel ~ ED_32L_left$leftMid, data = ED_32L, span=loessSpan)
loess_ED_32L_left <- predict(smooth_ED_32L_left)
smooth_ED_32L_right <- loess(ED_32L_right$Vel ~ ED_32L_right$rightMid, data = ED_32L, span=loessSpan)
loess_ED_32L_right <- predict(smooth_ED_32L_right)
## kernal density 
ED_kd_32L <- ggplot(ED_32L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=ED_32L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
ED_32L_kelp <- ggplot(ED_32L_left, aes(x=MidPt, y=loess_ED_32L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
ED_32L_urchin <- ggplot(ED_32L_right, aes(x=MidPt, y=loess_ED_32L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
ED_32L_both <- ED_kd_32L + annotation_custom(ggplotGrob(ED_32L_kelp)) + annotation_custom(ggplotGrob(ED_32L_urchin))
print(ED_32L_both)


## subset based on direction
ED_39R_left <- filter(ED_39R, left < 0) 
ED_39R_right <- filter(ED_39R, right > 0)
## velocities
smooth_ED_39R_left <- loess(ED_39R_left$Vel ~ ED_39R_left$leftMid, data = ED_39R, span=loessSpan)
loess_ED_39R_left <- predict(smooth_ED_39R_left)
smooth_ED_39R_right <- loess(ED_39R_right$Vel ~ ED_39R_right$rightMid, data = ED_39R, span=loessSpan)
loess_ED_39R_right <- predict(smooth_ED_39R_right)
## kernal density 
ED_kd_39R <- ggplot(ED_39R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=ED_39R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
ED_39R_kelp <- ggplot(ED_39R_left, aes(x=MidPt, y=loess_ED_39R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
ED_39R_urchin <- ggplot(ED_39R_right, aes(x=MidPt, y=loess_ED_39R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
ED_39R_both <- ED_kd_39R + annotation_custom(ggplotGrob(ED_39R_kelp)) + annotation_custom(ggplotGrob(ED_39R_urchin))
print(ED_39R_both)


## subset based on direction
ED_45R_left <- filter(ED_45R, left < 0) 
ED_45R_right <- filter(ED_45R, right > 0)
## velocities
smooth_ED_45R_left <- loess(ED_45R_left$Vel ~ ED_45R_left$leftMid, data = ED_45R, span=loessSpan)
loess_ED_45R_left <- predict(smooth_ED_45R_left)
smooth_ED_45R_right <- loess(ED_45R_right$Vel ~ ED_45R_right$rightMid, data = ED_45R, span=loessSpan)
loess_ED_45R_right <- predict(smooth_ED_45R_right)
## kernal density 
ED1 <- grobTree(text_grob("East Dutch", x=.25, y=.75, hjust=0, size = 12, color = "black"))
ED_kd_45R <- ggplot(ED_45R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=ED_45R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes +
  annotation_custom(ED1)
## plot velocity lines
ED_45R_kelp <- ggplot(ED_45R_left, aes(x=MidPt, y=loess_ED_45R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
ED_45R_urchin <- ggplot(ED_45R_right, aes(x=MidPt, y=loess_ED_45R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
ED_45R_both <- ED_kd_45R + annotation_custom(ggplotGrob(ED_45R_kelp)) + annotation_custom(ggplotGrob(ED_45R_urchin))
print(ED_45R_both)
## END East Dutch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## West Dutch ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## subset based on direction
WD_10R_left <- filter(WD_10R, left < 0) 
WD_10R_right <- filter(WD_10R, right > 0)
## velocities
smooth_WD_10R_left <- loess(WD_10R_left$Vel ~ WD_10R_left$leftMid, data = WD_10R, span=loessSpan)
loess_WD_10R_left <- predict(smooth_WD_10R_left)
smooth_WD_10R_right <- loess(WD_10R_right$Vel ~ WD_10R_right$rightMid, data = WD_10R, span=loessSpan)
loess_WD_10R_right <- predict(smooth_WD_10R_right)
## kernal density 
WD_kd_10R <- ggplot(WD_10R, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WD_10R_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
WD_10R_kelp <- ggplot(WD_10R_left, aes(x=MidPt, y=loess_WD_10R_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WD_10R_urchin <- ggplot(WD_10R_right, aes(x=MidPt, y=loess_WD_10R_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WD_10R_both <- WD_kd_10R + annotation_custom(ggplotGrob(WD_10R_kelp)) + annotation_custom(ggplotGrob(WD_10R_urchin))
print(WD_10R_both)


## subset based on direction
WD_22L_left <- filter(WD_22L, left < 0) 
WD_22L_right <- filter(WD_22L, right > 0)
## velocities
smooth_WD_22L_left <- loess(WD_22L_left$Vel ~ WD_22L_left$leftMid, data = WD_22L, span=loessSpan)
loess_WD_22L_left <- predict(smooth_WD_22L_left)
smooth_WD_22L_right <- loess(WD_22L_right$Vel ~ WD_22L_right$rightMid, data = WD_22L, span=loessSpan)
loess_WD_22L_right <- predict(smooth_WD_22L_right)
## kernal density 
WD_kd_22L <- ggplot(WD_22L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WD_22L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
WD_22L_kelp <- ggplot(WD_22L_left, aes(x=MidPt, y=loess_WD_22L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WD_22L_urchin <- ggplot(WD_22L_right, aes(x=MidPt, y=loess_WD_22L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WD_22L_both <- WD_kd_22L + annotation_custom(ggplotGrob(WD_22L_kelp)) + annotation_custom(ggplotGrob(WD_22L_urchin))
print(WD_22L_both)


## subset based on direction
WD_32L_left <- filter(WD_32L, left < 0) 
WD_32L_right <- filter(WD_32L, right > 0)
## velocities
smooth_WD_32L_left <- loess(WD_32L_left$Vel ~ WD_32L_left$leftMid, data = WD_32L, span=loessSpan)
loess_WD_32L_left <- predict(smooth_WD_32L_left)
smooth_WD_32L_right <- loess(WD_32L_right$Vel ~ WD_32L_right$rightMid, data = WD_32L, span=loessSpan)
loess_WD_32L_right <- predict(smooth_WD_32L_right)
## kernal density 
WD_kd_32L <- ggplot(WD_32L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WD_32L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + 
  show.below + xlab("NMDS Axis-1: System State")
## plot velocity lines
WD_32L_kelp <- ggplot(WD_32L_left, aes(x=MidPt, y=loess_WD_32L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WD_32L_urchin <- ggplot(WD_32L_right, aes(x=MidPt, y=loess_WD_32L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WD_32L_both <- WD_kd_32L + annotation_custom(ggplotGrob(WD_32L_kelp)) + annotation_custom(ggplotGrob(WD_32L_urchin))
print(WD_32L_both)


## subset based on direction
WD_39L_left <- filter(WD_39L, left < 0) 
WD_39L_right <- filter(WD_39L, right > 0)
## velocities
smooth_WD_39L_left <- loess(WD_39L_left$Vel ~ WD_39L_left$leftMid, data = WD_39L, span=loessSpan)
loess_WD_39L_left <- predict(smooth_WD_39L_left)
smooth_WD_39L_right <- loess(WD_39L_right$Vel ~ WD_39L_right$rightMid, data = WD_39L, span=loessSpan)
loess_WD_39L_right <- predict(smooth_WD_39L_right)
## kernal density 
WD_kd_39L <- ggplot(WD_39L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WD_39L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes
## plot velocity lines
WD_39L_kelp <- ggplot(WD_39L_left, aes(x=MidPt, y=loess_WD_39L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WD_39L_urchin <- ggplot(WD_39L_right, aes(x=MidPt, y=loess_WD_39L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WD_39L_both <- WD_kd_39L + annotation_custom(ggplotGrob(WD_39L_kelp)) + annotation_custom(ggplotGrob(WD_39L_urchin))
print(WD_39L_both)


## subset based on direction
WD_45L_left <- filter(WD_45L, left < 0) 
WD_45L_right <- filter(WD_45L, right > 0)
## velocities
smooth_WD_45L_left <- loess(WD_45L_left$Vel ~ WD_45L_left$leftMid, data = WD_45L, span=loessSpan)
loess_WD_45L_left <- predict(smooth_WD_45L_left)
smooth_WD_45L_right <- loess(WD_45L_right$Vel ~ WD_45L_right$rightMid, data = WD_45L, span=loessSpan)
loess_WD_45L_right <- predict(smooth_WD_45L_right)
## kernal density 
WD1 <- grobTree(text_grob("West Dutch", x=.25, y=.75, hjust=0, size = 12, color = "black"))
WD_kd_45L <- ggplot(WD_45L, aes(NMDS1)) + geom_density(adjust = kernalAdjust, fill=fillcolor, color=WD_45L_col, lwd=lineWeight_kd) + 
  my.theme + scale_x_continuous(limits = c(xmin, xmax)) + scale_y_continuous(limits=c(ymin_kd, ymax_kd_day)) + no.axes +
  annotation_custom(WD1)
## plot velocity lines
WD_45L_kelp <- ggplot(WD_45L_left, aes(x=MidPt, y=loess_WD_45L_left)) + geom_line(color=loessLeft, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
WD_45L_urchin <- ggplot(WD_45L_right, aes(x=MidPt, y=loess_WD_45L_right)) + geom_line(color=loessRight, lwd=lineWeight_loess) + scale_x_continuous(limits = c(xmin, xmax)) +
  scale_y_continuous(limits = c(ymin_line, ymax_line), expand = c(0,0)) + theme_bw() + vel.theme
## plot both
WD_45L_both <- WD_kd_45L + annotation_custom(ggplotGrob(WD_45L_kelp)) + annotation_custom(ggplotGrob(WD_45L_urchin))
print(WD_45L_both)
## End West Dutch 





graphics.off()
windows(h=24,w=20, record = TRUE)


allfigs <- ggarrange(tag_facet(NF_10R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(NF_22R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(NF_32L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(NF_39R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(NF_45R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45R" ),
                 
                 
                 tag_facet(WEK_10R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(WEK_22R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(WEK_32L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(WEK_39R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(WEK_45L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 
                 tag_facet(WEU_10L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10L"),
                 tag_facet(WEU_22L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(WEU_32R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32R"),
                 tag_facet(WEU_39L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39L"),
                 tag_facet(WEU_45L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 
                 tag_facet(Day_10R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(Day_22L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(Day_22R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(Day_32L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(Day_39L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R" ),
                 
                 
                 tag_facet(ED_10R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(ED_22R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(ED_32L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(ED_39R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(ED_45R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45R" ),
                 
                 
                 tag_facet(WD_10R_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(WD_22L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(WD_32L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(WD_39L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39L"),
                 tag_facet(WD_45L_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 nrow=6)


print(allfigs)





































