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

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
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
x.min <- -1.6
x.max <- 1.6
x.scale <- scale_x_continuous(limits = c(x.min, x.max))
y.min.kd <- 0
y.max.kd <- 3
y.scale.kd <- scale_y_continuous(limits = c(y.min.kd, y.max.kd))

## kernal density options 
fill.col <- "#CDCDCD"
line.weight <- .75
kernal.adjust <- 1

## loess options
loess.span <- 0.75
loess.left <- "#008000"
loess.right <- "#660198"
loess.line.weight <- .75
y.min.line <- 0
y.max.line <- 0.01
y.scale.line <- scale_y_continuous(limits = c(y.min.line, y.max.line), expand = c(0,0))
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





## functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## calculate loess smoother
smooth.f <- function(y, x){ ## e.g. NF_10R_left$Vel, NF_10R_left$MidPt, OR, within apply.f(dat_L$Vel, dat_L$MidPt)
  smooth_L <- loess(y ~ x, span=loess.span)
  predict_L <- predict(smooth_L)
}


## plot kernal density 
plot.density <- function(dat, col){ ## e.g. NF_10R, NF_10R_col 
  ggplot(dat, aes(x=NMDS1)) + geom_density(adjust=kernal.adjust, fill=fill.col, color=col, lwd=line.weight) +
    my.theme + x.scale + y.scale.kd + no.axes
}


## plot loess-smoothed velocity 
plot.vel <- function(dat, mid, y, loess.col){ ## e.g. NF_10R_left, NF_10R_left$MidPt, smooth_L, loess.left 
  ggplot(dat, aes(mid, y)) + geom_line(color=loess.col, lwd=line.weight) + 
    theme_bw() + vel.theme + x.scale + y.scale.line
}


## apply previously defined functions to calculate loess-smoother velocity and plot 
apply.f <- function(transect, col){
dat_L <- filter(transect, left < 0)
dat_R <- filter(transect, right > 0)
smooth_L <- smooth.f(dat_L$Vel, dat_L$MidPt)
smooth_R <- smooth.f(dat_R$Vel, dat_R$MidPt)
den <- plot.density(transect, col)
vel_L <- plot.vel(dat_L, dat_L$MidPt, smooth_L, loess.left)
vel_R <- plot.vel(dat_R, dat_R$MidPt, smooth_R, loess.right)
out <- den + annotation_custom(ggplotGrob(vel_L)) + annotation_custom(ggplotGrob(vel_R))
}
## END functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## apply functions to create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p.NF_10R <- apply.f(NF_10R, NF_10R_col)
p.NF_22R <- apply.f(NF_22R, NF_22R_col)
p.NF_32L <- apply.f(NF_32L, NF_32L_col)
p.NF_39R <- apply.f(NF_39R, NF_39R_col)
p.NF_45R <- apply.f(NF_45R, NF_45R_col)

p.WEK_10R <- apply.f(WEK_10R, WEK_10R_col)
p.WEK_22R <- apply.f(WEK_22R, WEK_22R_col)
p.WEK_32L <- apply.f(WEK_32L, WEK_32L_col)
p.WEK_39R <- apply.f(WEK_39R, WEK_39R_col)
p.WEK_45L <- apply.f(WEK_45L, WEK_45L_col)

p.WEU_10L <- apply.f(WEU_10L, WEU_10L_col)
p.WEU_22L <- apply.f(WEU_22L, WEU_22L_col)
p.WEU_32R <- apply.f(WEU_32R, WEU_32R_col)
p.WEU_39L <- apply.f(WEU_39L, WEU_39L_col)
p.WEU_45L <- apply.f(WEU_45L, WEU_45L_col)

p.Day_10R <- apply.f(Day_10R, Day_10R_col)
p.Day_22L <- apply.f(Day_22L, Day_22L_col)
p.Day_22R <- apply.f(Day_22R, Day_22R_col)
p.Day_32L <- apply.f(Day_32L, Day_32L_col)
p.Day_39L <- apply.f(Day_39L, Day_39L_col)

p.ED_10R <- apply.f(ED_10R, ED_10R_col)
p.ED_22R <- apply.f(ED_22R, ED_22R_col)
p.ED_32L <- apply.f(ED_32L, ED_32L_col)
p.ED_39R <- apply.f(ED_39R, ED_39R_col)
p.ED_45R <- apply.f(ED_45R, ED_45R_col)

p.WD_10R <- apply.f(WD_10R, WD_10R_col)
p.WD_22L <- apply.f(WD_22L, WD_22L_col)
p.WD_32L <- apply.f(WD_32L, WD_32L_col)
p.WD_39L <- apply.f(WD_39L, WD_39L_col)
p.WD_45L <- apply.f(WD_45L, WD_45L_col)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=24,w=20, record = TRUE)


figS1 <- ggarrange(tag_facet(p.NF_10R + facet_wrap(~"NMDS1"), tag_pool = "10R"),
                     tag_facet(p.NF_22R + facet_wrap(~"NMDS1"), tag_pool = "22R"),
                     tag_facet(p.NF_32L + facet_wrap(~"NMDS1"), tag_pool = "32L"),
                     tag_facet(p.NF_39R + facet_wrap(~"NMDS1"), tag_pool = "39R"),
                     tag_facet(p.NF_45R + facet_wrap(~"NMDS1"), tag_pool = "45R" ),
                     tag_facet(p.WEK_10R + facet_wrap(~"NMDS1"), tag_pool = "10R"),
                     tag_facet(p.WEK_22R + facet_wrap(~"NMDS1"), tag_pool = "22R"),
                     tag_facet(p.WEK_32L + facet_wrap(~"NMDS1"), tag_pool = "32L"),
                     tag_facet(p.WEK_39R + facet_wrap(~"NMDS1"), tag_pool = "39R"),
                     tag_facet(p.WEK_45L + facet_wrap(~"NMDS1"), tag_pool = "45L" ),
                     tag_facet(p.WEU_10L + facet_wrap(~"NMDS1"), tag_pool = "10L"),
                     tag_facet(p.WEU_22L + facet_wrap(~"NMDS1"), tag_pool = "22L"),
                     tag_facet(p.WEU_32R + facet_wrap(~"NMDS1"), tag_pool = "32R"),
                     tag_facet(p.WEU_39L + facet_wrap(~"NMDS1"), tag_pool = "39L"),
                     tag_facet(p.WEU_45L + facet_wrap(~"NMDS1"), tag_pool = "45L" ),
                     tag_facet(p.Day_10R + facet_wrap(~"NMDS1"), tag_pool = "10R"),
                     tag_facet(p.Day_22L + facet_wrap(~"NMDS1"), tag_pool = "22L"),
                     tag_facet(p.Day_22R + facet_wrap(~"NMDS1"), tag_pool = "22R"),
                     tag_facet(p.Day_32L + facet_wrap(~"NMDS1"), tag_pool = "32L"),
                     tag_facet(p.Day_39L + facet_wrap(~"NMDS1"), tag_pool = "39R" ),
                     tag_facet(p.ED_10R + facet_wrap(~"NMDS1"), tag_pool = "10R"),
                     tag_facet(p.ED_22R + facet_wrap(~"NMDS1"), tag_pool = "22R"),
                     tag_facet(p.ED_32L + facet_wrap(~"NMDS1"), tag_pool = "32L"),
                     tag_facet(p.ED_39R + facet_wrap(~"NMDS1"), tag_pool = "39R"),
                     tag_facet(p.ED_45R + facet_wrap(~"NMDS1"), tag_pool = "45R" ),
                     tag_facet(p.WD_10R + facet_wrap(~"NMDS1"), tag_pool = "10R"),
                     tag_facet(p.WD_22L + facet_wrap(~"NMDS1"), tag_pool = "22L"),
                     tag_facet(p.WD_32L + facet_wrap(~"NMDS1"), tag_pool = "32L"),
                     tag_facet(p.WD_39L + facet_wrap(~"NMDS1"), tag_pool = "39L"),
                     tag_facet(p.WD_45L + facet_wrap(~"NMDS1"), tag_pool = "45L" ),
                     nrow=6)

print(figS1)
## END plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
