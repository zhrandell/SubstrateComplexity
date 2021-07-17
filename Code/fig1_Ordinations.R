



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

## rugosity data 
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
relief_dat <- read.csv("SubstrateRugosity.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## process rugosity values and add to ordination data sheet ~~~~~~~~~~~~~~~~~~~~
relief_dat$Rugosity <- relief_dat$RELIEF / 10


## calculate mean, sd, se for rugosity values 
new_rugosity <- ddply(relief_dat, c("SITE", "TRANSECT", "ID"), summarise, 
                      N = length(Rugosity),
                      mean_rug = mean(Rugosity),
                      sd_rug = sd(Rugosity),
                      se_rug = sd_rug / sqrt(N))


## filter out sandy cove
new_rug <- filter(new_rugosity, SITE %in% c("NavFac","West End Kelp","West End Urchin",
                                            "Daytona","East Dutch","West Dutch"))


## join rugosity values to ordination data 
newdat <- dplyr::inner_join(dat, new_rug, by="ID")

## fix names of original file columns (altered from table join)
names(newdat)[7]<-"SITE"
names(newdat)[10]<-"TRANSECT"

## delete redunant columns (added from table join)
dat <- newdat[, !(colnames(newdat) %in% c("SITE.y","TRANSECT.y"))]
## data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





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
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphing ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal_All <- c("#ab2221", "#be5e00", "#7b6800", "#4d6934", "#3386a2", "#803280")

graphics.off()
windows(h=6,w=6, record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Create Giant kelp ordination ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
K1 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.85, stroke = 0.0, aes(color=SITE, size=MacPyr)) + 
  coord_fixed() + scale_size_continuous(range = c(0.9,5.5)) + scale_color_manual(values=pal_All) + 
  theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.5, 1.1) + ylim(-1.1, 1.4) + ylab("NMDS Axis-2") +
  guides(color = FALSE, size = guide_legend("Giant Kelp", order = 1, fill="black")) + 
  theme(legend.position=c(1,1), legend.justification=c(.98,.98), legend.text=element_text(size=12), legend.title = element_text(size=12),
        legend.title.align = 0.5, legend.key.height = unit(.2, "cm"), legend.key=element_rect(fill = FALSE, color=FALSE), legend.background=element_rect(fill=FALSE, color=FALSE),
        axis.title.x = element_blank(), axis.title.y = element_text(size=13), axis.ticks = element_blank(), plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) 
print(K1)


## custom annotations for macroalgae window
Mp <- grobTree(text_grob("Mp", x=.34, y=.9, hjust=0, size = 10, color = "black", face = "bold"))
j.Mp <- grobTree(text_grob("j.Mp", x=.1, y=.86, hjust=0, size = 10, color = "black", face = "bold"))
Pc <- grobTree(text_grob("Pc", x=.05, y=.35, hjust=0, size = 10, color = "black", face = "bold"))
Ls <- grobTree(text_grob("Ls", x=.35, y=.5, hjust=0, size = 10, color = "black", face = "bold"))
Ea <- grobTree(text_grob("Ea", x=.55, y=.55, hjust=0, size = 10, color = "black", face = "bold"))
j.Ls <- grobTree(text_grob("j.Ls", x=.1, y=.6, hjust=0, size = 10, color = "black", face = "bold"))
So <- grobTree(text_grob("So", x=.39, y=.2, hjust=0, size = 10, color = "black", face = "bold"))


## Giant kelp window with macroalgae correlation coefficients
K5 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-.35, 0) + ylim(-.15, .075) +
  guides(color = FALSE, size = FALSE) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend= -0.17185619 , yend= -0.13675533), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.22435289 , yend= -0.07666191), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.25878253, yend= 0.04577237), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.26220114, yend=-0.01797976), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.14890798, yend=  -0.04760285), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.34203979, yend= -0.10737747), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.17897971, yend= 0.05760792), color="darkgreen", lwd=0.80,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  annotation_custom(Mp) +
  annotation_custom(j.Mp) +
  annotation_custom(Pc) +
  annotation_custom(Ls) +
  annotation_custom(Ea) +
  annotation_custom(j.Ls) +
  annotation_custom(So)
print(K5)


## combine Giant kelp ordination with correlation plot 
KM <- K1 + annotation_custom(ggplotGrob(K5), 
                             xmin= -1.6, xmax= -.6, 
                             ymin= -1.25, ymax= -0.5)

print(KM)
## END Giant Kelp ordination ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## Second pane with herbivores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
K2 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.85, stroke = 0.0, aes(color=SITE, size=StrPur)) + 
  coord_fixed() + scale_size_continuous(range = c(0.9,5.5)) + scale_color_manual(values=pal_All) + 
  theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.5, 1.1) + ylim(-1.1, 1.4) +
  guides(color = FALSE, size = guide_legend("Purple Urchin", order = 1, fill="black")) + 
  theme(legend.position=c(1,1), legend.justification=c(.98,.98), legend.text=element_text(size=12), 
        legend.title = element_text(size=12), legend.title.align = 0.5, legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE), legend.background=element_rect(fill=FALSE, color=FALSE)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) 
print(K2)

## extract legend from K4 and add to K2
K4 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.85, aes(color=SITE, size=mean_rug)) + 
  coord_fixed() + scale_size_continuous(range = c(1,4), breaks = c(10.5, 12.5, 15.0, 17.5, 20.0, 22.5)) + scale_color_manual(values=pal_All) + 
  theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.5, 1.1) + ylim(-1.1, 1.4) + ylab("Ordination Axis 2") +
  guides(color = guide_legend(override.aes = list(size=5), nrow=1),
         size = FALSE) + 
  theme(legend.text=element_text(size=12), legend.key=element_rect(fill = FALSE, color=FALSE), legend.background=element_rect(fill=FALSE, color=FALSE),
        legend.title=element_blank(), legend.position="bottom",
        axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x=0, y=0, xend=0.078041, yend=0.996950), color="black",lwd=.5, arrow = arrow(length = unit(0.02, "npc"))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 


## extract legend from existing plot --- K4 
leg2 <- gtable_filter(ggplot_gtable(ggplot_build(K4)), "guide-box")


## Legend custom placed on K2 
K2 <- K2 + annotation_custom(grob=leg2, xmin= -.25, xmax=0, ymin= 3.3, ymax=0)
print(K2)


## custom annotations for herbivore correlations 
Sp <- grobTree(text_grob("Sp", x=.55, y=.55, hjust=0, size = 10, color = "black", face="bold"))
Mf <- grobTree(text_grob("Mf", x=.5, y=.85, hjust=0, size = 10, color = "black", face="bold"))
Mu <- grobTree(text_grob("Mu", x=.8, y=.2, hjust=0, size = 10, color = "black", face="bold"))


## inset pane with herbivore correlation coefficients
K6 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) + coord_fixed() +
  theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(0, .25) + ylim(-.2, .1) +
  guides(color = FALSE, size = FALSE) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend=0.10254074, yend=0.06198384), color="#91219E",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= 0.11854649, yend= -0.03334409), color="#91219E",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= 0.20499612, yend= -0.18911673), color="#91219E",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) + 
  annotation_custom(Sp) +
  annotation_custom(Mf) +
  annotation_custom(Mu)
print(K6)


KU <- K2 + annotation_custom(ggplotGrob(K6), 
                             xmin= .45, xmax= 1.325, 
                             ymin= -1.2, ymax= -0.5)

print(KU)
## END second pane with herbivores ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## pane 3 -- sea star ordination w/ pt size as rugosity ~~~~~~~~~~~~~~~~~~~~~~~~
K3 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.85, stroke = 0.0, aes(color=SITE, size=mean_rug)) + 
  coord_fixed() + scale_size_continuous(range = c(0.9,4.5), breaks = c(1.1, 1.2, 1.5, 1.7, 2.0)) + scale_color_manual(values=pal_All) + 
  theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.5, 1.1) + ylim(-1.1, 1.4) +
  guides(color = FALSE, size = guide_legend("Rugosity", order = 1, fill="black")) + 
  theme(legend.position=c(1,1), legend.justification=c(.98,.98), legend.text=element_text(size=12), legend.title = element_text(size=12),
        legend.title.align = 0.5, legend.key.height = unit(.2, "cm"), legend.key=element_rect(fill = FALSE, color=FALSE), legend.background=element_rect(fill=FALSE, color=FALSE)) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  geom_segment(aes(x=0, y=0, xend=0.078041, yend=0.996950), color="black",lwd=1,
               arrow = arrow(length = unit(0.02, "npc"))) +
  theme(plot.margin = margin(r=.1, l=.1, b=2, t=20, unit = "pt")) 
print(K3)


## custom annotations 
Di <- grobTree(text_grob("Di", x=.3, y=.75, hjust=0, size = 10, color = "black", face = "bold"))
Pm <- grobTree(text_grob("Pm", x=0, y=.65, hjust=0, size = 10, color = "black", face = "bold"))
Pg <- grobTree(text_grob("Pg", x=.8, y=.575, hjust=0, size = 10, color = "black", face = "bold"))
Ph <- grobTree(text_grob("Ph", x=.6, y=.9, hjust=0, size = 10, color = "black", face = "bold"))


## plot inset pane with sea star correlations 
K7 <- ggplot(dat, aes(x = NMDS1, y = NMDS2)) +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-.10, .08) + ylim(-.005, .275) +
  guides(color = FALSE, size = FALSE) + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank()) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  geom_segment(aes(x=0, y=0, xend=-0.09241174, yend=0.16711635), color="#0072B2",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= -0.04561673, yend= 0.19262086), color="#0072B2",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= 0.05017742, yend= 0.13358665), color="#0072B2",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  geom_segment(aes(x=0, y=0, xend= 0.05867194, yend= 0.26269635), color="#0072B2",lwd=.75,
               arrow = arrow(angle = 35, length = unit(.3, "cm"))) +
  annotation_custom(Pm) +
  annotation_custom(Di) +
  annotation_custom(Pg) +
  annotation_custom(Ph) 
print(K7)


## combine ordination and pane 
KS <- K3 + annotation_custom(ggplotGrob(K7), 
                             xmin= .1, xmax= .8, 
                             ymin= 0.6, ymax= 1.5)

print(KS)
## END third pane with sea star correlations and rugosity overlay ~~~~~~~~~~~~~~





## combine all three ordinations fto create figure 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=6,w=18, record=TRUE)


fig1 <- ggarrange(tag_facet(KM + facet_wrap(~"NMDS1"), tag_pool = "a"),
                  tag_facet(KU + facet_wrap(~"NMDS1"), tag_pool="b"),
                  tag_facet(KS + facet_wrap(~"NMDS1"), tag_pool="c"),
                  nrow=1, bottom = "NMDS Axis-1", 
                  label.args=list(gp=grid::gpar(size=13)))

print(fig1)
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


        
