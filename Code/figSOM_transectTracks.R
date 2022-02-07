## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## FigS3 in SubstrateComplexity: individual transect trajectories through two 
## dimensional species space
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)
library(grid)
library(gridExtra)
library(ggpubr)
library(gtable)
library(egg)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)
## END initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"))


## filter by site
filter.site <- function(site){
  filter(dat, SITE %in% c(site))
}

NF <- filter.site("NavFac")
WEK <- filter.site("WestEnd Kelp")
WEU <- filter.site("WestEnd Urchin")
Day <- filter.site("Daytona")
ED <- filter.site("East Dutch")
WD <- filter.site("West Dutch")


## filter by transect
filter.transect <- function(site, transect){
  filter(site, TRANSECT %in% c(transect))
}

NF_10R <- filter.transect(NF, "10R")
NF_22R <- filter.transect(NF, "22R")
NF_32L <- filter.transect(NF, "32L")
NF_39R <- filter.transect(NF, "39R")
NF_45R <- filter.transect(NF, "45R")

WEK_10R <- filter.transect(WEK, "10R")
WEK_22R <- filter.transect(WEK, "22R")
WEK_32L <- filter.transect(WEK, "32L")
WEK_39R <- filter.transect(WEK, "39R")
WEK_45L <- filter.transect(WEK, "45L")

WEU_10L <- filter.transect(WEU, "10L")
WEU_22L <- filter.transect(WEU, "22L")
WEU_32R <- filter.transect(WEU, "32R")
WEU_39L <- filter.transect(WEU, "39L")
WEU_45L <- filter.transect(WEU, "45L")

Day_10R <- filter.transect(Day, "10R")
Day_22L <- filter.transect(Day, "22L")
Day_22R <- filter.transect(Day, "22R")
Day_32L <- filter.transect(Day, "32L")
Day_39L <- filter.transect(Day, "39L")

ED_10R <- filter.transect(ED, "10R")
ED_22R <- filter.transect(ED, "22R")
ED_32L<- filter.transect(ED, "32L")
ED_39R <- filter.transect(ED, "39R")
ED_45R <- filter.transect(ED, "45R")

WD_10R <- filter.transect(WD, "10R")
WD_22L <- filter.transect(WD, "22L")
WD_32L <- filter.transect(WD, "32L")
WD_39L <- filter.transect(WD, "39L")
WD_45L <- filter.transect(WD, "45L")


## custom cols to call within plotting function
NF.cols <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392") 
WEK.cols <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") 
WEU.cols <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") 
Day.cols <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") 
ED.cols <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1")
WD.cols <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7")


## create row titles (site names) to append onto figures
title.f <- function(x){
  grobTree(text_grob(x, x=0.5, y=0.06, size=13, color="black"))
  }

NF.title <- title.f("NavFac")
WEK.title <- title.f("WestEnd Kelp")
WEU.title <- title.f("WestEnd Urchin")
Day.title <- title.f("Daytona")
ED.title <- title.f("East Dutch")
WD.title <- title.f("West Dutch")
## END data wrangling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pt.type <- 20
alph <- 0.9
pt.size <- 0.25
line.size <- 0.1

x.scale <- xlim(-1.4, 1.0) 
y.scale <- ylim(-1.0, 1.4)
margin <- theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 


my.theme <- theme(panel.border = element_rect(color="gray50"), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(), 
                  legend.position="none")

no.titles <- theme(axis.title.x = element_blank(),
                   axis.title.y = element_blank())

x.title <- theme(axis.title.x = element_text(size=15),
                 axis.title.y = element_blank())

y.title <- theme(axis.title.x = element_blank(), 
                 axis.title.y = element_text(size=15))
## END graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphing functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.f <- function(dat, col, t){
  ggplot(dat, aes(NMDS1, NMDS2)) + 
    geom_point(pch=pt.type, alpha=alph, size=pt.size, color=col) + 
    geom_path(size=line.size, color="black") + 
    coord_fixed() + theme_bw() + my.theme + x.scale + y.scale + margin + 
    no.titles 
}


plot.xlab <- function(dat, col){
  ggplot(dat, aes(NMDS1, NMDS2)) + 
    geom_point(pch=pt.type, alpha=alph, size=pt.size, color=col) + 
    geom_path(size=line.size, color="black") + 
    coord_fixed() + theme_bw() + my.theme + x.scale + y.scale + margin + 
    x.title + xlab("NMDS Axis-1")
}


plot.ylab <- function(dat, col){
  ggplot(dat, aes(NMDS1, NMDS2)) + 
    geom_point(pch=pt.type, alpha=alph, size=pt.size, color=col) + 
    geom_path(size=line.size, color="black") + 
    coord_fixed() + theme_bw() + my.theme + x.scale + y.scale + margin + 
    y.title + ylab("NMDS Axis-2")
}
## END graphing functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NF10R <- plot.f(NF_10R, NF.cols[[1]])
NF22R <- plot.f(NF_22R, NF.cols[[2]]) 
NF32L <- plot.f(NF_32L, NF.cols[[3]]) 
NF39R <- plot.f(NF_39R, NF.cols[[4]]) 
NF45R <- plot.f(NF_45R, NF.cols[[5]]) 
NF45R <- NF45R + annotation_custom(NF.title)

WEK10R <- plot.f(WEK_10R, WEK.cols[[1]])
WEK22R <- plot.f(WEK_22R, WEK.cols[[2]]) 
WEK32L <- plot.f(WEK_32L, WEK.cols[[3]]) 
WEK39R <- plot.f(WEK_39R, WEK.cols[[4]]) 
WEK45L <- plot.f(WEK_45L, WEK.cols[[5]]) 
WEK45L <- WEK45L + annotation_custom(WEK.title)

WEU10L <- plot.ylab(WEU_10L, WEU.cols[[1]])
WEU22L <- plot.f(WEU_22L, WEU.cols[[2]]) 
WEU32R <- plot.f(WEU_32R, WEU.cols[[3]]) 
WEU39L <- plot.f(WEU_39L, WEU.cols[[4]]) 
WEU45L <- plot.f(WEU_45L, WEU.cols[[5]]) 
WEU45L <- WEU45L + annotation_custom(WEU.title)

Day10R <- plot.f(Day_10R, Day.cols[[1]])
Day22L <- plot.f(Day_22L, Day.cols[[2]]) 
Day22R <- plot.f(Day_22R, Day.cols[[3]]) 
Day32L <- plot.f(Day_32L, Day.cols[[4]]) 
Day39L <- plot.f(Day_39L, Day.cols[[5]]) 
Day39L <- Day39L + annotation_custom(Day.title)

ED10R <- plot.f(ED_10R, ED.cols[[1]])
ED22R <- plot.f(ED_22R, ED.cols[[2]]) 
ED32L <- plot.f(ED_32L, ED.cols[[3]]) 
ED39R <- plot.f(ED_39R, ED.cols[[4]]) 
ED45R <- plot.f(ED_45R, ED.cols[[5]]) 
ED45R <- ED45R + annotation_custom(ED.title)

WD10R <- plot.f(WD_10R, WD.cols[[1]])
WD22L <- plot.f(WD_22L, WD.cols[[2]]) 
WD32L <- plot.xlab(WD_32L, WD.cols[[3]]) 
WD39L <- plot.f(WD_39L, WD.cols[[4]]) 
WD45L <- plot.f(WD_45L, WD.cols[[5]]) 
WD45L <- WD45L + annotation_custom(WD.title)
## END plot creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=24,w=20, record = TRUE)

figS3 <- ggarrange(tag_facet(NF10R + facet_wrap(~"NMDS1"), tag_pool = "a"),
                   tag_facet(NF22R + facet_wrap(~"NMDS1"), tag_pool = "b"),
                   tag_facet(NF32L + facet_wrap(~"NMDS1"), tag_pool = "c"),
                   tag_facet(NF39R + facet_wrap(~"NMDS1"), tag_pool = "d"),
                   tag_facet(NF45R + facet_wrap(~"NMDS1"), tag_pool = "e" ),
                   tag_facet(WEK10R + facet_wrap(~"NMDS1"), tag_pool = "f"),
                   tag_facet(WEK22R + facet_wrap(~"NMDS1"), tag_pool = "g"),
                   tag_facet(WEK32L + facet_wrap(~"NMDS1"), tag_pool = "h"),
                   tag_facet(WEK39R + facet_wrap(~"NMDS1"), tag_pool = "i"),
                   tag_facet(WEK45L + facet_wrap(~"NMDS1"), tag_pool = "j"),
                   tag_facet(WEU10L + facet_wrap(~"NMDS1"), tag_pool = "k"),
                   tag_facet(WEU22L + facet_wrap(~"NMDS1"), tag_pool = "l"),
                   tag_facet(WEU32R + facet_wrap(~"NMDS1"), tag_pool = "m"),
                   tag_facet(WEU39L + facet_wrap(~"NMDS1"), tag_pool = "n"),
                   tag_facet(WEU45L + facet_wrap(~"NMDS1"), tag_pool = "o" ),
                   tag_facet(Day10R + facet_wrap(~"NMDS1"), tag_pool = "p"),
                   tag_facet(Day22L + facet_wrap(~"NMDS1"), tag_pool = "q"),
                   tag_facet(Day22R + facet_wrap(~"NMDS1"), tag_pool = "r"),
                   tag_facet(Day32L + facet_wrap(~"NMDS1"), tag_pool = "s"),
                   tag_facet(Day39L + facet_wrap(~"NMDS1"), tag_pool = "t" ),
                   tag_facet(ED10R + facet_wrap(~"NMDS1"), tag_pool = "u"),
                   tag_facet(ED22R + facet_wrap(~"NMDS1"), tag_pool = "v"),
                   tag_facet(ED32L + facet_wrap(~"NMDS1"), tag_pool = "w"),
                   tag_facet(ED39R + facet_wrap(~"NMDS1"), tag_pool = "x"),
                   tag_facet(ED45R + facet_wrap(~"NMDS1"), tag_pool = "y"),
                   tag_facet(WD10R + facet_wrap(~"NMDS1"), tag_pool = "x"),
                   tag_facet(WD22L + facet_wrap(~"NMDS1"), tag_pool = "a1"),
                   tag_facet(WD32L + facet_wrap(~"NMDS1"), tag_pool = "a2"),
                   tag_facet(WD39L + facet_wrap(~"NMDS1"), tag_pool = "a3"),
                   tag_facet(WD45L + facet_wrap(~"NMDS1"), tag_pool = "a4"),
                   nrow=6, ncol=5)
## END plot all ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
