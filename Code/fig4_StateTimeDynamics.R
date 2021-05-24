## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Script to produce Figure 4, State Dynamics over time ~~~~~~~~~~~~~~~~~~~~~ ##
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

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

graphics.off()
windows(h=10,w=6, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"), 
                   labels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                            "Daytona", "East Dutch", "West Dutch"))


## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plotting configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color palattes 
pal_NF <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392") #BE2625 used for tint and shade creation
pal_WEK <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") #FF6600 used for tint and shade creation
pal_WEU <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") #CDAD00 used for tint and shade creation
pal_Day <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") #608341 used for tint and shade creation
pal_ED <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1") #00688B used for tint and shade creation
pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD used for tint and shade creation
pal_SC <- c("#062620", "#106050", "#1a997f", "#36c5a9", "#79d9c5")


## custom annotations 
algae <- grobTree(text_grob("Algae", x=1, y=-5, hjust=0, size = 16, color = "#308014", face = "bold"))
mixed <- grobTree(text_grob("Mixed", x=0, y=-10, hjust=0, size = 16, color = "#FFA500", face = "bold"))
barren <- grobTree(text_grob("Barren", x=1, y=-15, hjust=0, size = 15, color = "#660198", face = "bold"))


## sequence to plot along axis 2 
yrz=seq(1981,2018)
## END Plotting configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Plot Navfac 
p1 <- ggplot(NF, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_NF) +
  geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  scale_y_reverse(limits=c(74,0),
                  breaks=seq(0,74, by=2),
                  labels = yrz) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_text(size=13),
        axis.ticks.x = element_blank()) +
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p1)





## Plot WestEnd Urchin 
p2 <- ggplot(WEU, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WEU) +
  geom_point() +
  geom_path() +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=74.5,label="WestEnd Urchin",size=4.75,vjust=2.2,hjust=.5)
print(p2)





## Plot WestEnd Kelp 
p3 <- ggplot(WEK, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WEK) +
  geom_point() +
  geom_path() +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=74.5,label="WestEnd Kelp",size=4.75,vjust=2.2,hjust=.5)
print(p3)





## Print Day 
p4 <- ggplot(Day, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_Day) +
  geom_point() +
  geom_path() +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  xlab("NMDS Axis-1: System State") +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=74.5,label="Daytona",size=4.75,vjust=2.2,hjust=.5)
print(p4)





## Print EastDutch 
p5 <- ggplot(ED, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_ED) +
  geom_point() +
  geom_path() +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=74.5,label="East Dutch",size=4.75,vjust=2.2,hjust=.5)
print(p5)





## Print WestDutch 
p6 <- ggplot(WD, aes(NMDS1,PERIOD,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WD) +
  geom_point() +
  geom_path() +
  scale_y_reverse() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  xlim(-1.25, 1.2) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank())+
  annotate("text", x=1,y=0,label="Barren",size=4.75,vjust=-.7,hjust=.6) +
  annotate("text", x=0,y=0,label="Mixed",size=4.75,vjust=-.7,hjust=.7) +
  annotate("text", x=-1,y=0,label="Algae",size=4.75,vjust=-.7,hjust=.5) +
  annotate("text", x=0,y=74.5,label="West Dutch",size=4.75,vjust=2.2,hjust=.5)
print(p6)




## final figure
windows(h=10,w=25, record=TRUE)

g1 <- grid.arrange(arrangeGrob(p1,p2,p3,p4,p5,p6, 
                               nrow=1, ncol=6),
                   bottom=textGrob("NMDS Axis-1: System State",
                                   gp=gpar(fontsize=14)))

print(g1)
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
