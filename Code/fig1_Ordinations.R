



## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(egg)
library(tidyverse)
library(grid)
library(ggpubr)
library(gtable)


## rugosity data 
setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

cc <- read.csv("spp_correlation_coefficients.csv", header = TRUE)
cc <- cc[,-1]
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## convert substrate linear relief -> rugosity 
dat$rugosity <- dat$mean_rug / 10


## reorder sites to plot desired order 
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin", 
                                      "Daytona", "East Dutch", "West Dutch")) 


## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal_All <- c("#ab2221", "#be5e00", "#b3a100", "#4d6934", "#3386a2", "#803280")
myCols <- scale_color_manual(values=pal_All)

myTheme <- theme(panel.border = element_rect(color="gray50"), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank()) 

paneTheme <- theme(axis.title.x = element_blank(), 
                   axis.title.y = element_blank(), 
                   axis.text.x = element_blank(), 
                   axis.text.y = element_blank(), 
                   axis.ticks = element_blank())

legendTheme <- theme(legend.position=c(1,1),
                     legend.justification=c(.98,.98), 
                     legend.text=element_text(size=13), 
                     legend.title = element_text(size=13), 
                     legend.title.align = 0.5, 
                     legend.key.height = unit(.2, "cm"),
                     legend.key=element_rect(fill = FALSE, color=FALSE), 
                     legend.background=element_rect(fill=FALSE, color=FALSE))

bottom <- theme(axis.title.x = element_blank(), 
                axis.title.y = element_blank(), 
                axis.text.y = element_blank(), 
                axis.ticks = element_blank(), 
                axis.text.x = element_text(size=12)) 
    
sizeMin <- 0.9
sizeMax <- 5.5 
ptSize <- scale_size_continuous(range = c(sizeMin, sizeMax))
xMin <- -1.5 
xMax <- 1.1
yMin <- -1.1
yMax <- 1.4

ptType = 20 
alph = 0.85
str = 0

hjust <- 0
paneText <- 10
paneCol <- "black"
paneFace <- "bold"

kelpCol <- "darkgreen"
urchinCol <- "#91219E"
starCol <- "#0072B2"
lineWeight <- 0.8
arrowLength <- 0.3
arrowUnit <- "cm"
arrowAngle <- 35

noLegs <- guides(color = FALSE, size = FALSE)
## END graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Create Giant kelp ordination ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=6,w=6, record=TRUE)


p1 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = ptType, alpha = alph, stroke = str, aes(color=SITE, size=MacPyr)) + coord_fixed() + 
  ptSize + myCols + theme_bw() + myTheme + xlim(xMin, xMax) + ylim(yMin, yMax) + ylab("NMDS Axis-2") +
  guides(color = FALSE, size = guide_legend("Giant Kelp", order = 1, fill="black")) + legendTheme + 
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), 
        axis.text = element_text(size=12), axis.ticks = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) +
  annotate("text", x=-1.5, y=1.4, label="A", size=8)

print(p1)

## custom annotations for macroalgae window
Mp <- grobTree(text_grob("Mp", x=.34, y=.9, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
j.Mp <- grobTree(text_grob("j.Mp", x=.1, y=.86, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Pc <- grobTree(text_grob("Pc", x=.05, y=.35, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Ls <- grobTree(text_grob("Ls", x=.35, y=.5, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Ea <- grobTree(text_grob("Ea", x=.55, y=.55, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
j.Ls <- grobTree(text_grob("j.Ls", x=.1, y=.6, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
So <- grobTree(text_grob("So", x=.39, y=.2, hjust=hjust, size = paneText, color = paneCol, face = paneFace))


## Giant kelp window with macroalgae correlation coefficients
p2 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() + theme_bw() + myTheme + noLegs + paneTheme +  
  xlim(-.35, 0) + ylim(-.15, .075) + 
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend= cc[14,2], yend= cc[14,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[10,2], yend= cc[10,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[13,2], yend= cc[13,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[9,2], yend= cc[9,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[12,2], yend= cc[12,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[11,2], yend= cc[11,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend= cc[8,2], yend= cc[8,3]), color=kelpCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  annotation_custom(Mp) + annotation_custom(j.Mp) + annotation_custom(Pc) + annotation_custom(Ls) +
  annotation_custom(Ea) + annotation_custom(j.Ls) + annotation_custom(So)


## combine Giant kelp ordination with correlation plot 
fig1a <- p1 + annotation_custom(ggplotGrob(p2), 
                             xmin= -1.6, xmax= -.6, 
                             ymin= -1.25, ymax= -0.5)

print(fig1a)
## END Giant Kelp ordination ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Second pane with herbivores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
p4 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = ptType, alpha = alph, stroke = str, aes(color=SITE, size=StrPur)) + 
  coord_fixed() + theme_bw() + ptSize + myCols + myTheme + xlim(xMin, xMax) + ylim(yMin, yMax) +
  guides(color = FALSE, size = guide_legend("Purple Urchin", order = 1, fill="black")) + 
  legendTheme + bottom + theme(plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) +
  annotate("text", x=-1.5, y=1.4, label="B", size=8)



## extract legend from p5 and add to p4
p5 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) + myCols + 
  geom_point(pch = ptType, alpha = alph, aes(color=SITE, size=mean_rug)) +  
  guides(color = guide_legend(override.aes = list(size=5), nrow=1), size = FALSE) + 
  theme(legend.text=element_text(size=14), legend.key=element_rect(fill = FALSE, color=FALSE), 
        legend.background=element_rect(fill=FALSE, color=FALSE), legend.title=element_blank(), 
        legend.position="bottom", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 


## extract legend from p5 and add to p4
leg1 <- gtable_filter(ggplot_gtable(ggplot_build(p5)), "guide-box")
p4 <- p4 + annotation_custom(grob=leg1, xmin= -.25, xmax=0, ymin= 3.2, ymax=0)


## custom annotations for herbivore correlations 
Sp <- grobTree(text_grob("Sp", x=.55, y=.55, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Mf <- grobTree(text_grob("Mf", x=.5, y=.85, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Mu <- grobTree(text_grob("Mu", x=.8, y=.2, hjust=hjust, size = paneText, color = paneCol, face = paneFace))


## inset pane with herbivore correlation coefficients
p5 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) + coord_fixed() +
  theme_bw() + myTheme + xlim(0, .25) + ylim(-.2, .1) +  noLegs + paneTheme + 
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend = cc[2,2], yend = cc[2,3]), color=urchinCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend = cc[3,2], yend = cc[3,3]), color=urchinCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend = cc[4,2], yend = cc[4,3]), color=urchinCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) + 
  annotation_custom(Sp) + annotation_custom(Mf) + annotation_custom(Mu)


## combine into final pane
fig1b <- p4 + annotation_custom(ggplotGrob(p5), 
                             xmin= .45, xmax= 1.325, 
                             ymin= -1.2, ymax= -0.5)

print(fig1b)
## END second pane with herbivores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## pane 3 -- sea star ordination w/ pt size as rugosity ~~~~~~~~~~~~~~~~~~~~~~~~
p6 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.85, stroke = 0.0, aes(color=SITE, size=rugosity)) + 
  coord_fixed() + scale_size_continuous(range = c(sizeMin, 4.5), breaks = c(1.01, 1.3, 1.6, 1.9, 2.2)) + 
  myCols +  theme_bw() + myTheme + bottom + xlim(xMin, xMax) + ylim(yMin, yMax) +
  guides(color = FALSE, size = guide_legend("Rugosity", order = 1, fill="black")) + legendTheme + 
  geom_segment(aes(x=0, y=0, xend=0.158, yend=0.987), color="black",lwd=1, arrow = arrow(length = unit(0.02, "npc"))) +
  theme(plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) +
  annotate("text", x=-1.5, y=1.4, label="C", size=8)



## custom annotations 
Di <- grobTree(text_grob("Di", x=.3, y=.75, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Pm <- grobTree(text_grob("Pm", x=0, y=.65, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Pg <- grobTree(text_grob("Pg", x=.8, y=.575, hjust=hjust, size = paneText, color = paneCol, face = paneFace))
Ph <- grobTree(text_grob("Ph", x=.6, y=.9, hjust=hjust, size = paneText, color = paneCol, face = paneFace))


## plot inset pane with sea star correlations 
p7 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-.10, .08) + ylim(-.005, .275) + noLegs + paneTheme + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend = cc[1,2], yend = cc[1,3]), color=starCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend = cc[5,2], yend = cc[5,3]), color=starCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend = cc[6,2], yend = cc[6,3]), color=starCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  geom_segment(aes(x=0, y=0, xend = cc[7,2], yend = cc[7,3]), color=starCol, lwd=lineWeight, arrow = arrow(angle = arrowAngle, length = unit(arrowLength, arrowUnit))) +
  annotation_custom(Pm) + annotation_custom(Di) + annotation_custom(Pg) + annotation_custom(Ph) 


## combine ordination and pane 
fig1c <- p6 + annotation_custom(ggplotGrob(p7), 
                                xmin= .18, xmax= .8, 
                                ymin= 0.6, ymax= 1.5)

print(fig1c)
## END third pane with sea star correlations and rugosity overlay ~~~~~~~~~~~~~~





## combine all three ordinations to create fig1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=6,w=18, record=TRUE)

fig1 <- ggarrange(fig1a, fig1b, fig1c, nrow=1)
annotate_figure(fig1, bottom=text_grob("NMDS Axis-1", size=14))
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        
