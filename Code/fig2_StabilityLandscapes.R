## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Figure 2 for SubstrateComplexity ms ~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~ ##
## Modified May 25th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## reorder sites for plotting 
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch")) 
                   

## subset by site
NF <- filter(dat, SITE %in% c("NavFac"))
WEK <- filter(dat, SITE %in% c("WestEnd Kelp"))
Day <- filter(dat, SITE %in% c("Daytona"))
ED <- filter(dat, SITE %in% c("East Dutch"))
WEU <- filter(dat, SITE %in% c("WestEnd Urchin"))
WD <- filter(dat, SITE %in% c("West Dutch"))


## subset data to plot individual transects 
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


## subset data to plot grayed out points in row 1 transect tracks 
ab_WEK_39R <- filter(dat, ID %in% c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                    "2_10L", "2_22", "2_32R", "2_39L", "2_45L",
                                    "4_10R", "4_22L", "4_32L", "4_39L", "4_45L",
                                    "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                    "6_10R", "6_22L", "6_22R", "6_32L", "6_39L"))


ab_NF_39R <- filter(dat, ID %in% c("2_10L", "2_22", "2_32R", "2_39L", "2_45L",
                                   "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                   "4_10R", "4_22L", "4_32L", "4_39L", "4_45L",
                                   "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                   "6_10R", "6_22L", "6_22R", "6_32L", "6_39L"))


ab_Day_22R <- filter(dat, ID %in% c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                    "2_10L", "2_22", "2_32R", "2_39L", "2_45L",
                                    "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                    "4_10R", "4_22L", "4_32L", "4_39L", "4_45L",
                                    "5_10R", "5_22R", "5_32L", "5_39R", "5_45R"))


ab_ED_22R <- filter(dat, ID %in% c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                   "2_10L", "2_22", "2_32R", "2_39L", "2_45L",
                                   "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                   "4_10R", "4_22L", "4_32L", "4_39L", "4_45L",
                                   "6_10R", "6_22L", "6_22R", "6_32L", "6_39L"))


ab_WD_45L <- filter(dat, ID %in% c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                   "2_10L", "2_22", "2_32R", "2_39L", "2_45L",
                                   "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                   "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                   "6_10R", "6_22L", "6_22R", "6_32L", "6_39L"))


ab_WEU_32R <- filter(dat, ID %in% c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                    "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                    "4_10R", "4_22L", "4_32L", "4_39L", "4_45L",
                                    "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                    "6_10R", "6_22L", "6_22R", "6_32L", "6_39L"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color palattes
pal_NF <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392") #BE2625 used for tint and shade creation
pal_WEK <- c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66") #FF6600 used for tint and shade creation
pal_WEU <- c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680") #CDAD00 used for tint and shade creation
pal_Day <- c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d") #608341 used for tint and shade creation
pal_ED <- c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1") #00688B used for tint and shade creation
pal_WD <- c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD used for tint and shade creation
#pal_SC <- c("#062620", "#106050", "#1a997f", "#36c5a9", "#79d9c5")


## custom annotations 
NF_grob <- grobTree(text_grob("NavFac 39R", x=.45, y=.05, hjust=0, size = 11, color = "black"))
WEK_grob <- grobTree(text_grob("WestEnd Kelp 39R", x=.45, y=.05, hjust=0, size = 11, color = "black"))
WEU_grob <- grobTree(text_grob("WestEnd Urchin 32R", x=.45, y=.05, hjust=0, size = 11, color = "black"))
Day_grob <- grobTree(text_grob("Daytona 22R", x=.45, y=.05, hjust=0, size = 11, color = "black"))
ED_grob <- grobTree(text_grob("East Dutch 22R", x=.45, y=.05, hjust=0, size = 10, color = "black"))
WD_grob <- grobTree(text_grob("West Dutch 45L", x=.45, y=.05, hjust=0, size = 10, color = "black"))


graphics.off()
windows(w=4,h=4,record=TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## Plot Transect Tracks row 2 of figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom annotation 
Nf <- grobTree(text_grob("39R", x=.85, y=.95, hjust=0, size = 12, color = "black"))

NF_tr <- ggplot(ab_NF_39R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = NF_39R, pch = 20, alpha = 0.90, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_NF) + 
  geom_path(data = NF_39R, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) + ylab("NMDS Axis-2") +
  theme(axis.title.x = element_blank(), axis.title.y = element_text(size=14), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(),
        legend.position = "none") + theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Nf)

print(NF_tr)



Weu <- grobTree(text_grob("32R", x=.85, y=.95, hjust=0, size = 12, color = "black"))

WEU_tr <- ggplot(ab_WEU_32R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = WEU_32R, pch = 20, alpha = 0.90, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_WEU) + 
  geom_path(data = WEU_32R, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  xlab("System State") +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(-1.5, 1.2), position = "top") +
  ylim(-1.0, 1.4) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Weu)

print(WEU_tr)



Wek <- grobTree(text_grob("39R", x=.85, y=.95, hjust=0, size = 12, color = "black"))

WEK_tr <- ggplot(ab_WEK_39R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = WEK_39R, pch = 20, alpha = 0.90, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_WEK) + 
  geom_path(data = WEK_39R, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Wek)

print(WEK_tr)



Dy <- grobTree(text_grob("22R", x=.85, y=.95, hjust=0, size = 12, color = "black"))

Day_tr <- ggplot(ab_Day_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = Day_22R, pch = 20, alpha = 0.80, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_Day) + 
  geom_path(data = Day_22R, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Dy)

print(Day_tr)




Ed <- grobTree(text_grob("22R", x=.85, y=.95, hjust=0, size = 12, color = "black"))

ED_tr <- ggplot(ab_ED_22R, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = ED_22R, pch = 20, alpha = 0.90, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_ED) + 
  geom_path(data = ED_22R, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Ed)

print(ED_tr)




Wd <- grobTree(text_grob("45L", x=.85, y=.95, hjust=0, size = 12, color = "black"))

WD_tr <- ggplot(ab_WD_45L, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color="black", pch = 20, alpha = 0.05, fill="gray", size=.25) +
  geom_point(data = WD_45L, pch = 20, alpha = 0.90, aes(color=TRANSECT), size = .25) + 
  scale_color_manual(values=pal_WD) + 
  geom_path(data = WD_45L, aes(x=NMDS1, y=NMDS2), size=0.2) +
  coord_fixed() +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlim(-1.4, 1.0) + ylim(-1.0, 1.4) +
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), 
        legend.position = "none", plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt")) +
  annotation_custom(Wd)

print(WD_tr)
## END transect track plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~







## Plot kernal densities, velocities, and combine both ~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac Kernal densities 
NF_kd <- ggplot(NF, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 1, position="stack", color=NA) +
  scale_fill_manual(values=pal_NF) +
  scale_x_continuous(limits = c(-1.5, 1.2), position = "top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0)) +
  theme_bw() +
  ylab("Negative Potential") +
  xlab("NavFac") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14),
        axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(),
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"),
        legend.position="none") 
print(NF_kd)

## loess smoother through velocities
smooth_NF <- loess(NF$Vel ~ NF$MidPt, data = NF, span=.75)
loess_NF <- predict(smooth_NF)

## plot velocities and smoother 
Fin_NF <- ggplot(NF, aes(x=MidPt, y=loess_NF)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_NF)

## print both
NF_both <- NF_kd +
  annotation_custom(ggplotGrob(Fin_NF))
print(NF_both)





## West End Kelp kernal densities  
WEK_kd <- ggplot(WEK, aes(NMDS1, group=TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = .8, position = "stack", color=NA) +
  scale_fill_manual(values=pal_WEK) +
  scale_x_continuous(limits = c(-1.5, 1.2), position = "top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  xlab("West End Kelp") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"), legend.position="none") 
print(WEK_kd)

## loess smoother through velocities
smooth_WEK <- loess(WEK$Vel ~ WEK$MidPt, data = WEK, span=.75)
loess_WEK <- predict(smooth_WEK)

## plot velocities and smoother 
Fin_WEK <- ggplot(WEK, aes(x=MidPt, y=loess_WEK)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_WEK)

## print both
WEK_both <- WEK_kd +
  annotation_custom(ggplotGrob(Fin_WEK))
print(WEK_both)





## WestEnd Urchin kernal density
WEU_kd <- ggplot(WEU, aes(NMDS1, group=TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 0.8, position="stack", color=NA) +
  scale_fill_manual(values=pal_WEU) +
  scale_x_continuous(limits = c(-1.5, 1.2), position = "top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0)) +
  theme_bw() +
  xlab("West End Urchin") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"), legend.position="none") 
print(WEU_kd)

## loess smoother through velocities
smooth_WEU <- loess(WEU$Vel ~ WEU$MidPt, data = WEU, span=.75)
loess_WEU <- predict(smooth_WEU)

## plot velocities and smoother 
Fin_WEU <- ggplot(WEU, aes(x=MidPt, y=loess_WEU)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_WEU)

## print both
WEU_both <- WEU_kd +
  annotation_custom(ggplotGrob(Fin_WEU))
print(WEU_both)






## Daytona kernal density 
Day_kd <- ggplot(Day, aes(NMDS1, group=TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 1.0, position="stack", color=NA) +
  scale_fill_manual(values=pal_Day) +
  scale_x_continuous(limits = c(-1.5, 1.2), position="top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0)) +
  theme_bw() +
  xlab("Daytona") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"), legend.position="none") 
print(Day_kd)

## loess smoother through velocities
smooth_Day <- loess(Day$Vel ~ Day$MidPt, data = Day, span=.75)
loess_Day <- predict(smooth_Day)

## plot velocities and smoother 
Fin_Day <- ggplot(Day, aes(x=MidPt, y=loess_Day)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_Day)

## print both
Day_both <- Day_kd +
  annotation_custom(ggplotGrob(Fin_Day))
print(Day_both)

## East Dutch 
ED_kd <- ggplot(ED, aes(NMDS1, group=TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 1.0, position="stack", color=NA) +
  scale_fill_manual(values=pal_ED) +
  scale_x_continuous(limits = c(-1.5, 1.2), position="top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0)) +
  theme_bw() +
  xlab("East Dutch") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"), legend.position="none") 
print(ED_kd)

## loess smoother through velocities 
smooth_ED <- loess(ED$Vel ~ ED$MidPt, data = ED, span=.75)
loess_ED <- predict(smooth_ED)

## plot velocities and smoother 
Fin_ED <- ggplot(ED, aes(x=MidPt, y=loess_ED)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0)) +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA), 
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y.left = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_ED)

## print both 
ED_both <- ED_kd +
  annotation_custom(ggplotGrob(Fin_ED))
print(ED_both)







## West Dutch KERNAL DENSITY
WD_kd <- ggplot(WD, aes(NMDS1, group=TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 1.0, position="stack", color=NA) +
  scale_fill_manual(values=pal_WD) +
  scale_x_continuous(limits = c(-1.5, 1.2), position="top") +
  scale_y_continuous(limits = c(0,9.15), expand = c(0,0), position="right") +
  theme_bw() +
  xlab("West Dutch") +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  ylab("Velocity") +
  theme(axis.title.x = element_text(size=14), axis.title.y = element_text(size=14), axis.text.x = element_blank(), axis.text.y = element_blank(), axis.ticks = element_blank(), plot.title = element_blank(), 
        plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"), legend.position="none") 
print(WD_kd)

## loess smoother through velocities 
smooth_WD <- loess(WD$Vel ~ WD$MidPt, data = WD, span=.75)
loess_WD <- predict(smooth_WD)

## plot velocities and smoother 
Fin_WD <- ggplot(WD, aes(x=MidPt, y=loess_WD)) +
  geom_line(color = "black", lwd = 1.5) +
  geom_line(color = "white", lwd = .75) +
  theme(legend.position="none") +
  scale_x_continuous(limits = c(-1.5, 1.2)) +
  scale_y_continuous(limits = c(0,.01), expand = c(0,0), position="right") +
  theme_bw() +
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "transparent",colour = NA), plot.background = element_rect(fill = "transparent",colour = NA),
        plot.title = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank(), axis.text.y = element_blank(), axis.text.x = element_blank(), axis.ticks = element_blank()) 
print(Fin_WD)

## print both 
WD_both <- WD_kd +
  annotation_custom(ggplotGrob(Fin_WD))
print(WD_both)





## check velocities and kernal densities 
graphics.off()
windows(h=8,w=48, record = TRUE)


tr3 <- ggarrange(tag_facet(NF_both + facet_wrap(~"NMDS1"), tag_pool = "a"),
                 tag_facet(WEK_both + facet_wrap(~"NMDS1"), tag_pool = "b"),
                 tag_facet(WEU_both + facet_wrap(~"NMDS1"), tag_pool = "c"),
                 tag_facet(Day_both + facet_wrap(~"NMDS1"), tag_pool = "d"),
                 tag_facet(ED_both + facet_wrap(~"NMDS1"), tag_pool = "e" ),
                 tag_facet(WD_both + facet_wrap(~"NMDS1"), tag_pool = "f"),
                 nrow=1)
print(tr3)
## END kernal density and velocity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## LOAD substrate rugosity data for row 3 of figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~
library(ggplot2)
library(dplyr)
library(vegan)
library(MASS)
library(fitdistrplus)
library(egg)
library(plyr)
library(forcats)

setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
relief_dat <- read.csv("SubstrateRugosity.csv", header = TRUE)

graphics.off()
windows(h=6,w=6, record=TRUE)


## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## convert relief --> rugosity 
relief_dat$Rugosity <- relief_dat$RELIEF / 10


## calculate mean, sd, se 
new_rugosity <- ddply(relief_dat, c("SITE", "TRANSECT"), summarise, 
                      N = length(Rugosity),
                      mean_rug = mean(Rugosity),
                      sd_rug = sd(Rugosity),
                      se_rug = sd_rug / sqrt(N))


## filter data by site to plot 
NF_rug <- filter(new_rugosity, SITE %in% c("NavFac"))
WEK_rug <- filter(new_rugosity, SITE %in% c("West End Kelp"))
Day_rug <- filter(new_rugosity, SITE %in% c("Daytona"))
ED_rug <- filter(new_rugosity, SITE %in% c("East Dutch"))
WEU_rug <- filter(new_rugosity, SITE %in% c("West End Urchin"))
WD_rug <- filter(new_rugosity, SITE %in% c("West Dutch"))
#SC_rug <- filter(new_rugosity, SITE %in% c("Sandy Cove"))


## reorder transects (for plotting) by ascending values of rugosity 
NF_rug$TRANSECT <- factor(NF_rug$TRANSECT, levels = NF_rug$TRANSECT[order(NF_rug$mean_rug)])
WEK_rug$TRANSECT <- factor(WEK_rug$TRANSECT, levels = WEK_rug$TRANSECT[order(WEK_rug$mean_rug)])
WEU_rug$TRANSECT <- factor(WEU_rug$TRANSECT, levels = WEU_rug$TRANSECT[order(WEU_rug$mean_rug)])
Day_rug$TRANSECT <- factor(Day_rug$TRANSECT, levels = Day_rug$TRANSECT[order(Day_rug$mean_rug)])
ED_rug$TRANSECT <- factor(ED_rug$TRANSECT, levels = ED_rug$TRANSECT[order(ED_rug$mean_rug)])
WD_rug$TRANSECT <- factor(WD_rug$TRANSECT, levels = WD_rug$TRANSECT[order(WD_rug$mean_rug)])
#SC_rug$TRANSECT <- factor(SC_rug$TRANSECT, levels = SC_rug$TRANSECT[order(SC_rug$mean_rug)])
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom color palettes 
up_pal_SC <- c("#062620", "#36c5a9", "#106050", "#1a997f", "#79d9c5")
up_pal_NF <- c("#9f1009", "#5b0905", "#e3170d", "#e9453d", "#f18b86") #E3170D used for tint and shade creation
up_pal_WEK <- c("#b34700", "#662900", "#ff6600", "#ff8533", "#ffb380") #FF6600 used for tint and shade creation
up_pal_WEU <- c("#665c00", "#ffe600", "#b3a100", "#332e00", "#fff380") #FFE600 used for tint and shade creation
up_pal_Day <- c("#91c591", "#0a2a0a", "#228b22", "#4ea24e", "#145314")
up_pal_ED <- c("#99c3d1", "#003446", "#3386a2", "#00536f", "#66a4b9") #00688B used for tint and shade creation
up_pal_WD <- c("#d6adeb", "#b870dc", "#4d1967", "#7a28a4", "#a347d2") #9932CD used for tint and shade creation
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plotting ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac
NF1 <- ggplot(NF_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_NF) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank(), 
        axis.text.x = element_text(size=13, color = "black"), 
        axis.title.y = element_text(size = 14), 
        axis.text.y = element_text(size=13), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") +
  labs(x="NavFac") + 
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(NF1)


## West End Kelp 
WEK1 <- ggplot(WEK_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_WEK) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") +
  labs(x="West End Kelp") +
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(WEK1)


## West End Urchin 
WEU1 <- ggplot(WEU_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_WEU) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 14),
        axis.text.x = element_text(size=13, color = "black"),
        axis.title.y = element_blank(),
        axis.text.y = element_blank(), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") +
  labs(x="Transect") +
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(WEU1)


## Daytona
Day1 <- ggplot(Day_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_Day) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") + 
  labs(x="Daytona") +
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(Day1)


## East Dutch 
ED1 <- ggplot(ED_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_ED) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") + 
  labs(x="East Dutch") +
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(ED1)


## West Dutch 
WD1 <- ggplot(WD_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
  scale_fill_manual(values=up_pal_WD) + 
  theme_bw() +
  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.title.x = element_blank(),
        axis.text.x = element_text(size=13, color = "black"), 
        axis.title.y = element_blank(), 
        axis.text.y = element_blank(), 
        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
        legend.position="none") + 
  labs(x="West Dutch") +
  labs(y="Rugosity") +
  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
print(WD1)


## Plot Sandy Cove
#SC1 <- ggplot(SC_rug, aes(SITE, mean_rug)) +
#  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=.4) +
#  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = 3.5, shape=21, color="black") +
#  scale_fill_manual(values=up_pal_SC) + 
#  theme_bw() +
#  theme(panel.border = element_rect(color="gray50"), panel.grid.major = element_blank(),
#        panel.grid.minor = element_blank(), 
#        axis.title.x = element_blank(),
#        axis.text.x = element_text(size=13, color = "black"), 
#        axis.title.y = element_blank(), 
#        axis.text.y = element_blank(), 
#        plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
#        legend.position="none") + 
#  labs(x="Sandy Cove") +
#  labs(y="Rugosity") +
#  scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 
#print(SC1)
## End Site Plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## Aggregate all plots into single figure ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=4,w=24, record=TRUE)


n11 <- ggarrange(NF1, 
                 WEK1,
                 WEU1,
                 Day1,
                 ED1,
                 WD1,
                 nrow=1)

print(n11)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~






## combine all rows into final figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(w=24,h=12,record=TRUE)



tr3 <- ggarrange(tag_facet(NF_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "a"),
                 tag_facet(WEK_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "b"),
                 tag_facet(WEU_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "c"),
                 tag_facet(Day_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "d"),
                 tag_facet(ED_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "e" ),
                 tag_facet(WD_both +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "f"),
                 tag_facet(NF_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "g"),
                 tag_facet(WEK_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "h"),
                 tag_facet(WEU_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "i"),
                 tag_facet(Day_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "j"),
                 tag_facet(ED_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "k"),
                 tag_facet(WD_tr +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "l"),
                 tag_facet(NF1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "m"),
                 tag_facet(WEK1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "n"),
                 tag_facet(WEU1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "o"),
                 tag_facet(Day1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "p"),
                 tag_facet(ED1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "q"),
                 tag_facet(WD1 +
                             facet_wrap(~"NMDS1"),
                           tag_pool = "r"),
                 nrow=3, ncol=6)


print(tr3)
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##


                 
                 














