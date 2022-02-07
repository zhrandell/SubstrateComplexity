## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Figure 2 for SubstrateComplexity ms ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Modified Jan 31 2022; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)
library(MASS)
library(fitdistrplus)
library(gridExtra)
library(ggpubr)
library(egg)
library(gtable)
library(grid)
library(magick)
library(here)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
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


## subset data to the specific transect visualized in row 2 
NF_39R <- filter(NF, TRANSECT %in% c("39R")) 
WEK_39R <- filter(WEK, TRANSECT %in% c("39R")) 
WEU_32R <- filter(WEU, TRANSECT %in% c("32R")) 
WD_45L <- filter(WD, TRANSECT %in% c("45L"))
ED_22R <- filter(ED, TRANSECT %in% c("22R"))
Day_22R <- filter(Day, TRANSECT %in% c("22R"))


## subset data to plot grayed out points in row 1 transect tracks 
ab_NF <- filter(dat, SITE != c("NavFac"))
ab_WEK <- filter(dat, SITE != c("WestEnd Kelp"))
ab_WEU <- filter(dat, SITE != c("WestEnd Urchin"))
ab_Day <- filter(dat, SITE != c("Daytona"))
ab_ED <- filter(dat, SITE != c("East Dutch"))
ab_WD <- filter(dat, SITE != c("West Dutch"))
## END data subset ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## window size
graphics.off()
windows(w=4,h=4,record=TRUE)

## custom color palattes
NF.cols <- scale_color_manual(values=c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392"))
WEK.cols <- scale_color_manual(values=c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66"))
WEU.cols <- scale_color_manual(values=c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680"))
Day.cols <- scale_color_manual(values=c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d"))
ED.cols <- scale_color_manual(values=c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1"))
WD.cols <- scale_color_manual(values=c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7"))

## custom annotations 
NF.grob <- grobTree(text_grob("39R", x=.85, y=.95, hjust=0, size = 15, color = "black"))
WEU.grob <- grobTree(text_grob("32R", x=.85, y=.95, hjust=0, size = 15, color = "black"))
WEK.grob <- grobTree(text_grob("39R", x=.85, y=.95, hjust=0, size = 15, color = "black"))
Day.grob <- grobTree(text_grob("22R", x=.85, y=.95, hjust=0, size = 15, color = "black"))
ED.grob <- grobTree(text_grob("22R", x=.85, y=.95, hjust=0, size = 15, color = "black"))
WD.grob <- grobTree(text_grob("45L", x=.85, y=.95, hjust=0, size = 15, color = "black"))

## ggplot custom themes
bg.theme <- theme(panel.border = element_rect(color="gray50"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank())

my.theme <- theme(panel.border = element_rect(color="gray50"), 
                  panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(),
                  axis.text.x = element_blank(), 
                  axis.text.y = element_blank(), 
                  axis.ticks = element_blank(), 
                  legend.position="none") 

x.title <- theme(axis.title.x = element_text(size=15))
y.title <- theme(axis.title.y = element_text(size=15))
no.x.title <- theme(axis.title.x = element_blank())
no.y.title <- theme(axis.title.y = element_blank())
no.titles <- theme(axis.title = element_blank())
margins <- theme(plot.margin = margin(r=.1, l=.1, b=.1, t=.1, unit = "pt"))

## background (bg) pts 
bg.pt.type <- 20
bg.pt.alpha <- 0.05
bg.pt.fill <- "gray"
bg.pt.col <- "black"
bg.pt.size <- 0.25

## primary pts and line 
pt.type <- 20
pt.alpha <- 0.90
pt.size <- 0.25
path.size <- 0.20

## axis lims 
x.min <- -1.5
x.max <- 1.2
y.min <- -1.0
y.max <- 1.4 
x.axis <- xlim(x.min, x.max)
y.axis <- ylim(y.min, y.max)


## create sub-figure labels 
pane.lab.R2 <- function(x){
  annotate("text", x = -1.44, y=1.325, label=x, size = 8)
}

NF_tr.lab <- pane.lab.R2("G")
WEK_tr.lab <- pane.lab.R2("H")
WEU_tr.lab <- pane.lab.R2("I")
Day_tr.lab <- pane.lab.R2("J")
ED_tr.lab <- pane.lab.R2("K")
WD_tr.lab <- pane.lab.R2("L")
## END custom params for Fig2 row 2, transect tracks ~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot Transect Tracks row 2 of figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
NF_tr <- ggplot(ab_NF, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = NF_39R, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = NF_39R, aes(x=NMDS1, y=NMDS2), size = path.size) +
  NF.cols + coord_fixed() + theme_bw() + bg.theme + x.axis + y.axis + margins + my.theme + y.title + no.x.title + 
  ylab("NMDS Axis-2") + annotation_custom(NF.grob) + NF_tr.lab

WEU_tr <- ggplot(ab_WEU, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = WEU_32R, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = WEU_32R, aes(x=NMDS1, y=NMDS2), size = path.size) +
  scale_x_continuous(limits = c(-1.5, 1.2), position = "top") +
  WEU.cols + coord_fixed() + theme_bw() + bg.theme + y.axis + margins + my.theme + no.y.title + x.title +  
  xlab("System State") + annotation_custom(WEU.grob) + WEU_tr.lab

WEK_tr <- ggplot(ab_WEK, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = WEK_39R, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = WEK_39R, aes(x=NMDS1, y=NMDS2), size = path.size) +
  WEK.cols + coord_fixed() + theme_bw() + bg.theme + x.axis + y.axis + margins + my.theme + no.y.title + x.title + 
  xlab("NMDS Axis-1") + annotation_custom(WEK.grob) + WEK_tr.lab

Day_tr <- ggplot(ab_Day, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = Day_22R, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = Day_22R, aes(x=NMDS1, y=NMDS2), size = path.size) +
  Day.cols + coord_fixed() + theme_bw() + bg.theme + x.axis + y.axis + margins + my.theme + no.titles + 
  annotation_custom(Day.grob) + Day_tr.lab

ED_tr <- ggplot(ab_ED, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = ED_22R, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = ED_22R, aes(x=NMDS1, y=NMDS2), size = path.size) +
  ED.cols + coord_fixed() + theme_bw() + bg.theme + x.axis + y.axis + margins + my.theme + no.titles + 
  annotation_custom(ED.grob) + ED_tr.lab

WD_tr <- ggplot(ab_WD, aes(x = NMDS1, y = NMDS2)) +
  geom_point(color=bg.pt.col, pch = bg.pt.type, alpha = bg.pt.alpha, fill=bg.pt.fill, size=bg.pt.size) +
  geom_point(data = WD_45L, pch = pt.type, alpha = pt.alpha, aes(color=TRANSECT), size = pt.size) + 
  geom_path(data = WD_45L, aes(x=NMDS1, y=NMDS2), size = path.size) +
  WD.cols + coord_fixed() + theme_bw() + bg.theme + x.axis + y.axis + margins + my.theme + no.titles + 
  annotation_custom(WD.grob) + WD_tr.lab
## END transect track plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot kernal densities, velocities, and combine both ~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac Kernal densities 
NF.fill <- scale_fill_manual(values=c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392"))
WEK.fill <- scale_fill_manual(values=c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66"))
WEU.fill <- scale_fill_manual(values=c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680"))
Day.fill <- scale_fill_manual(values=c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d"))
ED.fill <- scale_fill_manual(values=c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1"))
WD.fill <- scale_fill_manual(values=c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7"))

loess.theme <- theme(panel.border = element_blank(), 
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank(), 
                     panel.background = element_rect(fill = "transparent",colour = NA), 
                     plot.background = element_rect(fill = "transparent",colour = NA), 
                     plot.title = element_blank(), 
                     axis.title.y = element_blank(), 
                     axis.title.x = element_blank(), 
                     axis.text.y.left = element_blank(), 
                     axis.text.x = element_blank(), 
                     axis.ticks = element_blank(), 
                     legend.position = "none")

den.adjust <- 1 
den.y.lim <- 9.15
den.y.axis <- scale_y_continuous(limits = c(0, den.y.lim), expand = c(0,0))
den.x.axis <- scale_x_continuous(limits = c(x.min, x.max), position = "top") 
den.margins <- theme(plot.margin = margin(r=.1, l=.1, b=.1, t=5, unit = "pt"))
loess.span <- 0.75

black.line <- geom_line(color = "black", lwd = 1.5)
white.line <- geom_line(color = "white", lwd = .75)
loess.ylim <- scale_y_continuous(limits = c(0,.01), expand = c(0,0))

## create sub-figure labels 
pane.lab.R1 <- function(x){
  annotate("text", x = -1.44, y=8.555, label=x, size = 8)
}

NF.R1 <- pane.lab.R1("A")
WEK.R1 <- pane.lab.R1("B")
WEU.R1 <- pane.lab.R1("C")
Day.R1 <- pane.lab.R1("D")
ED.R1 <- pane.lab.R1("E")
WD.R1 <- pane.lab.R1("F")
## END graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## velocity of community shift~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
apply.loess <- function(x, y, data){
  out1 <- loess(x ~ y, data, span = loess.span)
  out2 <- predict(out1)
}

## calculate velocities 
loess_NF <- apply.loess(NF$Vel, NF$MidPt, NF) 
loess_WEK <- apply.loess(WEK$Vel, WEK$MidPt, WEK) 
loess_WEU <- apply.loess(WEU$Vel, WEU$MidPt, WEU) 
loess_Day <- apply.loess(Day$Vel, Day$MidPt, Day) 
loess_ED <- apply.loess(ED$Vel, ED$MidPt, ED) 
loess_WD <- apply.loess(WD$Vel, WD$MidPt, WD) 
## END velocities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot kernal densities and velocity of community shifts (Fig2 row1) ~~~~~~~~~~
## NavFac kernal densities 
NF_kd <- ggplot(NF, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = den.adjust, position="stack", color=NA) + ylab("Negative Potential") + xlab("NavFac") +
  theme_bw() + NF.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + y.title + NF.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_NF <- ggplot(NF, aes(x=MidPt, y=loess_NF)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
NF_both <- NF_kd + annotation_custom(ggplotGrob(vel_NF))


## WestEnd Kelp kernal densities 
WEK_kd <- ggplot(WEK, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 0.8, position="stack", color=NA) + xlab("West End Kelp") +
  theme_bw() + WEK.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + no.y.title + WEK.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_WEK <- ggplot(WEK, aes(x=MidPt, y=loess_WEK)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
WEK_both <- WEK_kd + annotation_custom(ggplotGrob(vel_WEK))


## WestEnd Urchin kernal densities 
WEU_kd <- ggplot(WEU, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = 0.8, position="stack", color=NA) + xlab("West End Urchin") +
  theme_bw() + WEU.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + no.y.title + WEU.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_WEU <- ggplot(WEU, aes(x=MidPt, y=loess_WEU)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
WEU_both <- WEU_kd + annotation_custom(ggplotGrob(vel_WEU))


## Daytona kernal densities 
Day_kd <- ggplot(Day, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = den.adjust, position="stack", color=NA) + xlab("Daytona") +
  theme_bw() + Day.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + no.y.title + Day.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_Day <- ggplot(Day, aes(x=MidPt, y=loess_Day)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
Day_both <- Day_kd + annotation_custom(ggplotGrob(vel_Day))


## East Dutch kernal densities 
ED_kd <- ggplot(ED, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = den.adjust, position="stack", color=NA) + xlab("East Dutch") +
  theme_bw() + ED.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + no.y.title + ED.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_ED <- ggplot(ED, aes(x=MidPt, y=loess_ED)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
ED_both <- ED_kd + annotation_custom(ggplotGrob(vel_ED))


## West Dutch kernal densities 
WD_kd <- ggplot(WD, aes(NMDS1, group = TRANSECT, fill=TRANSECT)) +
  geom_density(adjust = den.adjust, position="stack", color=NA) + xlab("West Dutch") +
  theme_bw() + WD.fill + den.x.axis + den.y.axis + my.theme + den.margins + x.title + no.y.title + WD.R1
## plot velocities and combine kernal densities and velocities into single figure 
vel_WD <- ggplot(WD, aes(x=MidPt, y=loess_WD)) + 
  black.line + white.line + theme_bw() + x.axis + loess.ylim + loess.theme 
WD_both <- WD_kd + annotation_custom(ggplotGrob(vel_WD))
## END kernal density and velocity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## load substrate rugosity data for fig2 row 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
library(plyr)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
relief_dat <- read.csv("SubstrateRugosity.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## convert relief --> rugosity 
relief_dat$Rugosity <- relief_dat$RELIEF / 10


## calculate mean, sd, se 
rugosity <- ddply(relief_dat, c("SITE", "TRANSECT", "ID"), summarise, 
                  N = length(Rugosity),
                  mean = mean(Rugosity),
                  sd = sd(Rugosity),
                  se = sd / sqrt(N))


## remove Sandy Cove
rugosity <- filter(rugosity, SITE != c("Sandy Cove"))


## filter data by site to plot 
filter.site <- function(x){
  filter(rugosity, SITE %in% x)
}


NF_rug <- filter.site(c("NavFac"))
WEK_rug <- filter.site(c("West End Kelp"))
WEU_rug <- filter.site(c("West End Urchin"))
Day_rug <- filter.site(c("Daytona"))
ED_rug <- filter.site(c("East Dutch"))
WD_rug <- filter.site(c("West Dutch"))


## reorder transects (for plotting) by ascending values of rugosity 
NF_rug$TRANSECT <- factor(NF_rug$TRANSECT, levels = NF_rug$TRANSECT[order(NF_rug$mean)])
WEK_rug$TRANSECT <- factor(WEK_rug$TRANSECT, levels = WEK_rug$TRANSECT[order(WEK_rug$mean)])
WEU_rug$TRANSECT <- factor(WEU_rug$TRANSECT, levels = WEU_rug$TRANSECT[order(WEU_rug$mean)])
Day_rug$TRANSECT <- factor(Day_rug$TRANSECT, levels = Day_rug$TRANSECT[order(Day_rug$mean)])
ED_rug$TRANSECT <- factor(ED_rug$TRANSECT, levels = ED_rug$TRANSECT[order(ED_rug$mean)])
WD_rug$TRANSECT <- factor(WD_rug$TRANSECT, levels = WD_rug$TRANSECT[order(WD_rug$mean)])
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rug.theme <- theme(panel.border = element_rect(color="gray50"), 
                   panel.grid.major = element_blank(), 
                   panel.grid.minor = element_blank(),
                   axis.title.x = element_text(size=15),
                   axis.title.y = element_text(size=15),
                   axis.text.x = element_text(size=14), 
                   axis.text.y = element_text(size=14), 
                   plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
                   legend.position="none")

rug.y.scale <- scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 

no.y.title <- theme(axis.title.y = element_blank(), 
                    axis.text.y = element_blank())

no.x.title <- theme(axis.title.x = element_blank())

pt.size <- 3.5
pt.shape <- 21
pt.col <- "black"
sd.width <- 0.4


pane.lab.R3 <- function(x){
  annotate("text", x = .745, y=2.35, label=x, size = 8)
}

NF.R3 <- pane.lab.R3("M")
WEK.R3 <- pane.lab.R3("N")
WEU.R3 <- pane.lab.R3("O")
Day.R3 <- pane.lab.R3("P")
ED.R3 <- pane.lab.R3("Q")
WD.R3 <- pane.lab.R3("R")
## END custom params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## NavFac
p.NF <- ggplot(NF_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + NF.fill + rug.y.scale + labs(y="Rugosity") + no.x.title + NF.R3  

## WestEnd Kelp 
p.WEK <- ggplot(WEK_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + WEK.fill + rug.y.scale + no.y.title + no.x.title + WEK.R3

## WestEnd Urchin
p.WEU <- ggplot(WEU_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + WEU.fill + rug.y.scale + no.y.title + xlab("Transect") + WEU.R3

## Daytona
p.Day <- ggplot(Day_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + Day.fill + rug.y.scale + no.y.title + no.x.title + Day.R3

## East Dutch
p.ED <- ggplot(ED_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + ED.fill + rug.y.scale + no.y.title + no.x.title + ED.R3

## West Dutch
p.WD <- ggplot(WD_rug, aes(SITE, mean)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean-se, ymax=mean+se), width=sd.width) +
  geom_point(aes(x=TRANSECT, y=mean, fill=TRANSECT), size = pt.size, shape=pt.shape, color=pt.col) +
  theme_bw() + rug.theme + WD.fill + rug.y.scale + no.y.title + no.x.title + WD.R3
## END plots for fig2 row 3 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## combine all rows into final figure 2 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(w=24,h=12,record=TRUE)

fig2 <- ggarrange(NF_both, WEK_both, WEU_both, Day_both, ED_both, WD_both, 
                  NF_tr, WEK_tr, WEU_tr, Day_tr, ED_tr, WD_tr, 
                  p.NF, p.WEK, p.WEU, p.Day, p.ED, p.WD, 
                  nrow=3, ncol=6)
## END of plotting all subfigs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
