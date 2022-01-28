## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Script to produce Fig3 in "Kelp-forest dynamics controlled by substrate complexity" 
## updated January 28th 2022; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(egg)
library(tidyverse)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$SITE <- factor(dat$SITE, levels=c("NavFac", 
                                      "WestEnd Kelp", 
                                      "WestEnd Urchin",  
                                      "Daytona", 
                                      "East Dutch", 
                                      "West Dutch"))


## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))


## sequence to plot along axis 2 
years=seq(1980,2018)
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal_NF <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392") 
pal_WEK <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") 
pal_WEU <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") 
pal_Day <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") 
pal_ED <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1") 
pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") 
pal_SC <- c("#062620", "#106050", "#1a997f", "#36c5a9", "#79d9c5")


myTheme <- theme(panel.border = element_rect(color="gray50"), 
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank()) 

textTheme <- theme(axis.title.y = element_blank(),
                   axis.title.x = element_blank(),
                   axis.text.x = element_blank(),
                   axis.text.y = element_blank(),
                   axis.ticks = element_blank())

yAxis <- theme(axis.title.y = element_blank(),
               axis.title.x = element_blank(),
               axis.text.x = element_blank(),
               axis.text.y = element_text(size=16),
               axis.ticks.x = element_blank()) 


xMin <- -1.32
xMax <- 1.2 
yMin <- 0
yMax <- 76
noLeg <- guides(color=FALSE)
textSize <- 5.75


barren <- annotate("text", x=1, y=0, label= "Barren", size = textSize, vjust= -.7, hjust= .6) 
mixed <- annotate("text", x=0, y=0, label= "Mixed", size = textSize, vjust= -.7, hjust= .7) 
algae <- annotate("text", x=-1, y=0, label= "Algae", size = textSize, vjust= -.7, hjust= .5) 
## END plotting configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  




## Create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Navfac 
p1 <- ggplot(NF, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_NF) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + yAxis + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(a)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="NavFac", size = textSize, vjust=2.2, hjust=.5) 


## WestEnd Kelp 
p2 <- ggplot(WEK, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_WEK) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + textTheme + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(b)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="WestEnd Kelp", size = textSize, vjust=2.2, hjust=.5) 


## WestEnd Urchin 
p3 <- ggplot(WEU, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_WEU) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + textTheme + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(c)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="WestEnd Urchin", size = textSize, vjust=2.2, hjust=.5) 


## Daytona 
p4 <- ggplot(Day, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_Day) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + textTheme + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(d)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="Daytona", size = textSize, vjust=2.2, hjust=.5) 


## East Dutch 
p5 <- ggplot(ED, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_ED) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + textTheme + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(e)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="East Dutch", size = textSize, vjust=2.2, hjust=.5) 


## West Dutch 
p6 <- ggplot(WD, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale_color_manual(values=pal_WD) + xlim(xMin, xMax) +
  scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years) +
  theme_bw() + myTheme + noLeg + textTheme + barren + mixed + algae + 
  annotate("text", x = -1.2, y=75.5, label="(f)", size = textSize, vjust=2.2, hjust=.9) +
  annotate("text", x = 0, y=75.5, label="West Dutch", size = textSize, vjust=2.2, hjust=.5) 


## final figure
graphics.off()
windows(h=10,w=25, record=TRUE)

fig3 <- ggarrange(p1, p2, p3, p4, p5, p6, nrow=1)
annotate_figure(fig3, bottom=text_grob("NMDS Axis-1: System State", size=15))
## END figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
