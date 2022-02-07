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
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"))


## subset by site
filter.site <- function(x){
  filter(dat, SITE %in% c(x))
}

NF <- filter.site("NavFac")
WEK <- filter.site("WestEnd Kelp")
WEU <- filter.site("WestEnd Urchin")
Day <- filter.site("Daytona")
ED <- filter.site("East Dutch")
WD <- filter.site("West Dutch")


## sequence to plot along axis 2 
years=seq(1980,2018)
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NF.col <- scale_color_manual(values=c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392")) 
WEK.col <- scale_color_manual(values=c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66")) 
WEU.col <- scale_color_manual(values=c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680")) 
Day.col <- scale_color_manual(values=c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d")) 
ED.col <- scale_color_manual(values=c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1")) 
WD.col <- scale_color_manual(values=c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7")) 


## create sub-figure labels 
pane.lab <- function(x){
  annotate("text", x = -1.25, y=2, label=x, size = 8)
  }

NF.lab <- pane.lab("A")
WEK.lab <- pane.lab("B")
WEU.lab <- pane.lab("C")
Day.lab <- pane.lab("D")
ED.lab <- pane.lab("E")
WD.lab <- pane.lab("F")


## create site labels
site.lab <- function(site){
  annotate("text", x = 0, y=75, label=site, size = 6, vjust=2.2, hjust=.5)
  }

NF.site <- site.lab("NavFac")
WEK.site <- site.lab("WestEnd Kelp")
WEU.site <- site.lab("WestEnd Urchin")
Day.site <- site.lab("Daytona")
ED.site <- site.lab("East Dutch")
WD.site <- site.lab("West Dutch")


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
text.size <- 6
scale.y <- scale_y_reverse(limits=c(yMax, yMin), breaks=seq(yMin, yMax, by=2), labels = years)
scale.x <- xlim(xMin, xMax)

barren <- annotate("text", x=1, y=0, label= "Barren", size = text.size, vjust= -.7, hjust= .6) 
mixed <- annotate("text", x=0, y=0, label= "Mixed", size = text.size, vjust= -.7, hjust= .7) 
algae <- annotate("text", x=-1, y=0, label= "Algae", size = text.size, vjust= -.7, hjust= .5) 
## END plotting configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Navfac 
p1 <- ggplot(NF, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + yAxis + barren + mixed + algae + NF.lab + NF.site + NF.col

## WestEnd Kelp 
p2 <- ggplot(WEK, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + textTheme + barren + mixed + algae + WEK.lab + WEK.site + WEK.col

## WestEnd Urchin 
p3 <- ggplot(WEU, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + textTheme + barren + mixed + algae + WEU.lab + WEU.site + WEU.col

## Daytona 
p4 <- ggplot(Day, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + textTheme + barren + mixed + algae + Day.lab + Day.site + Day.col

## East Dutch 
p5 <- ggplot(ED, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + textTheme + barren + mixed + algae + ED.lab + ED.site + ED.col

## West Dutch 
p6 <- ggplot(WD, aes(NMDS1,PERIOD,color=TRANSECT)) +
  geom_point() + geom_path() + scale.y + scale.x + theme_bw() + myTheme + 
  noLeg + textTheme + barren + mixed + algae + WD.lab + WD.site + WD.col


## final figure
graphics.off()
windows(h=10,w=25, record=TRUE)

fig3 <- ggarrange(p1, p2, p3, p4, p5, p6, nrow=1)
annotate_figure(fig3, bottom=text_grob("NMDS Axis-1: System State", size=17))
## END figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
