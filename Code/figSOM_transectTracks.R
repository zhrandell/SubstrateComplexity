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





pal_NF <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392") #BE2625 used for tint and shade creation
pal_WEK <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") #FF6600 used for tint and shade creation
pal_WEU <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") #CDAD00 used for tint and shade creation
pal_Day <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") #608341 used for tint and shade creation
pal_ED <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1") #00688B used for tint and shade creation
pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD used for tint and shade creation



graphics.off()
windows(h=4,w=4, record = TRUE)



NF10R_tr <- ggplot(NF_10R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#390b0b") + 
  geom_path(size=0.1, color="black") +
  scale_color_manual(values=pal_NF) + coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(NF10R_tr)



NF22R_tr <- ggplot(NF_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90, size = .25, color="#721716") + 
  geom_path(size=0.1) +
  scale_color_manual(values=pal_NF) + coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(NF22R_tr)



NF32L_tr <- ggplot(NF_32L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#ab2221") + 
  geom_path(size=0.1) +
  scale_color_manual(values=pal_NF) + coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(NF32L_tr)



NF39R_tr <- ggplot(NF_39R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#cb5151") + 
  geom_path(size=0.1) +
  scale_color_manual(values=pal_NF) + coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(NF39R_tr)



NF5 <- grobTree(text_grob("NavFac", x=.35, y=.10, hjust=0, size = 12, color = "black"))
NF45R_tr <- ggplot(NF_45R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#df9392") + 
  geom_path(size=0.1) +
  scale_color_manual(values=pal_NF) + coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(),axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.2, unit = "pt")) +
  annotation_custom(NF5)

print(NF45R_tr)



## WEK 
pal_WEK <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") #FF6600 used for tint and shade creation

WEK10R_tr <- ggplot(WEK_10R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#301800") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEK10R_tr)



WEK22R_tr <- ggplot(WEK_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#773b00") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEK22R_tr)



WEK32L_tr <- ggplot(WEK_32L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#be5e00") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEK32L_tr)



WEK39R_tr <- ggplot(WEK_39R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#f0841a") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEK39R_tr)



WEK5 <- grobTree(text_grob("WestEnd Kelp", x=.15, y=.10, hjust=0, size = 12, color = "black"))
WEK45L_tr <- ggplot(WEK_45L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#f5ad66") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(WEK5)

print(WEK45L_tr)



pal_WEU <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") #CDAD00 used for tint and shade creation



WEU10L_tr <- ggplot(WEU_10L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#141100") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_text(size=13), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEU10L_tr)


WEU22L_tr <- ggplot(WEU_22L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#3d3400") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEU22L_tr)


WEU32R_tr <- ggplot(WEU_32R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#7b6800") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEU32R_tr)


WEU39L_tr <- ggplot(WEU_39L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#b99c00") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WEU39L_tr)


WEU5 <- grobTree(text_grob("WestEnd Urchin", x=.10, y=.10, hjust=0, size = 12, color = "black"))
WEU45L_tr <- ggplot(WEU_45L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#e6d680") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(WEU5)

print(WEU45L_tr)


pal_Day <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") #608341 used for tint and shade creation



Day10R_tr <- ggplot(Day_10R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#131a0d") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(Day10R_tr)


Day22L_tr <- ggplot(Day_22L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#304221") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(Day22L_tr)


Day22R_tr <- ggplot(Day_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#4d6934") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(Day22R_tr)

Day32L_tr <- ggplot(Day_32L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#708f54") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(Day32L_tr)


Day5 <- grobTree(text_grob("Daytona", x=.35, y=.10, hjust=0, size = 12, color = "black"))
Day39L_tr <- ggplot(Day_39L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#a0b58d") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Day5)

print(Day39L_tr)


pal_ED <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1") #00688B used for tint and shade creation
pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD used for tint and shade creation


ED10R_tr <- ggplot(ED_10R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#003446") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(ED10R_tr)



ED22R_tr <- ggplot(ED_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#00536f") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(ED22R_tr)



ED32L_tr <- ggplot(ED_32L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#3386a2") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(ED32L_tr)



ED39R_tr <- ggplot(ED_39R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#66a4b9") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(ED39R_tr)



ED5 <- grobTree(text_grob("East Dutch", x=.25, y=.10, hjust=0, size = 12, color = "black"))
ED45R_tr <- ggplot(ED_45R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#99c3d1") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(ED5)

print(ED45R_tr)




pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD used for tint and shade creation


WD10R_tr <- ggplot(WD_10R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#2b112b") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WD10R_tr)


WD22L_tr <- ggplot(WD_22L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#552255") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WD22L_tr)


WD32L_tr <- ggplot(WD_32L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#803280") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + xlab("NMDS Axis-1") +
  theme(axis.title.x = element_text(size=13), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WD32L_tr)


WD39L_tr <- ggplot(WD_39L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#a560a5") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), axis.text.x = element_blank(),axis.text.y = element_blank(),axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) 

print(WD39L_tr)


WD5 <- grobTree(text_grob("West Dutch", x=.25, y=.10, hjust=0, size = 12, color = "black"))
WD45L_tr <- ggplot(WD_45L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(pch = 20, alpha = 0.90,size = .25, color="#a560a5") + 
  geom_path(size=0.1, color="black") +
  coord_fixed() + theme_bw() + theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") + xlab("NMDS Axis-1") +
  theme(axis.title.x = element_blank(), 
        axis.title.y = element_blank(), 
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank()) +
  theme(legend.position=c(1.025,1.175), legend.justification=c(1,1), legend.text=element_text(size=7),legend.box = "horizontal", legend.title = element_blank(),
        legend.spacing.y = unit(.5, "cm"),
        legend.spacing.x = unit(0.0001, "cm"),
        legend.key.height = unit(.2, "cm"),
        legend.key=element_rect(fill = FALSE, color=FALSE),legend.background=element_rect(fill=FALSE, color=FALSE)) +
  guides(color = guide_legend(override.aes = list(size=2.5))) +
  theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(WD5)

print(WD45L_tr)











graphics.off()
windows(h=24,w=20, record = TRUE)


tr3 <- ggarrange(tag_facet(NF10R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(NF22R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(NF32L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(NF39R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(NF45R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45R" ),
                 
                 
                 tag_facet(WEK10R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(WEK22R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(WEK32L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(WEK39R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(WEK45L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 
                 tag_facet(WEU10L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10L"),
                 tag_facet(WEU22L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(WEU32R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32R"),
                 tag_facet(WEU39L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39L"),
                 tag_facet(WEU45L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 
                 tag_facet(Day10R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(Day22L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(Day22R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(Day32L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(Day39L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R" ),
                 
                 
                 tag_facet(ED10R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(ED22R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22R"),
                 tag_facet(ED32L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(ED39R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39R"),
                 tag_facet(ED45R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45R" ),
                 
                 
                 tag_facet(WD10R_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "10R"),
                 tag_facet(WD22L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "22L"),
                 tag_facet(WD32L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "32L"),
                 tag_facet(WD39L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "39L"),
                 tag_facet(WD45L_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "45L" ),
                 
                 nrow=6)


print(tr3)

