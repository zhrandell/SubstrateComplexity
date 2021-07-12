## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Script to produce SOM sheephead figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
dat <- read.csv("CAsheephead.csv", header = TRUE)

graphics.off()
windows(h=10,w=6, record=TRUE)

my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16), 
                 legend.text = element_text(size=14), 
                 legend.title = element_blank(), 
                 legend.position = c(0.85, 0.85))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data for kernal densities and eCDF ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$Site <- factor(dat$Station, 
                   levels=c("1","2","3","4","5","6"),
                   labels=c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))

## once more to reorder the site labels (for plotting)
dat$Site <- factor(dat$Site, levels=c("NavFac","West End Urchin","West End Kelp","Daytona","East Dutch","West Dutch"))

dat <- filter(dat, Site %in% c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))

dat$log <- log10(dat$SumAdultCount+1)


NF <- filter(dat, Site %in% c("NavFac"))
WEK <- filter(dat, Site %in% c("West End Kelp"))
WEU <- filter(dat, Site %in% c("West End Urchin"))
ED <- filter(dat, Site %in% c("East Dutch"))
WD <- filter(dat, Site %in% c("West Dutch"))
Day <- filter(dat, Site %in% c("Daytona"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot kernal densities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal_All <- c("#ab2221", "#be5e00", "#7b6800", "#4d6934", "#3386a2", "#803280")

graphics.off()
windows(w=8,h=5,record=TRUE)

p3 <- ggplot(dat, aes(log, fill=Site)) +
  geom_density(position="identity", color="black", alpha=0.3) +  my.theme +
  scale_fill_manual(values=pal_All) +
  xlab("total sheephead per site (per 0.1 hectare)") +
  guides(fill=guide_legend(order=1)) +
  scale_x_continuous(labels=c("0","3.16","10","31.6","100","316.2")) +
  theme(legend.title = element_blank())
print(p3)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## eCDFs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~
## calcuate empirical cumulative density function and extract and sort values
NF_ecdf <- data.frame(x=unique(NF$log), 
                      y=ecdf(NF$log)(unique(NF$log))*length(NF$log))
WEU_ecdf <- data.frame(x=unique(WEU$log), 
                       y=ecdf(WEU$log)(unique(WEU$log))*length(WEU$log))
WEK_ecdf <- data.frame(x=unique(WEK$log), 
                       y=ecdf(WEK$log)(unique(WEK$log))*length(WEK$log))
Day_ecdf <- data.frame(x=unique(Day$log), 
                       y=ecdf(Day$log)(unique(Day$log))*length(Day$log))
ED_ecdf <- data.frame(x=unique(ED$log), 
                      y=ecdf(ED$log)(unique(ED$log))*length(ED$log))
WD_ecdf <- data.frame(x=unique(WD$log), 
                      y=ecdf(WD$log)(unique(WD$log))*length(WD$log))


## rescale extracted cdf values to 0-1 scale
NF_ecdf$y <- scale(NF_ecdf$y, center=min(NF_ecdf$y), scale=diff(range(NF_ecdf$y)))
WEU_ecdf$y <- scale(WEU_ecdf$y, center=min(WEU_ecdf$y), scale=diff(range(WEU_ecdf$y)))
WEK_ecdf$y <- scale(WEK_ecdf$y, center=min(WEK_ecdf$y), scale=diff(range(WEK_ecdf$y)))
Day_ecdf$y <- scale(Day_ecdf$y, center=min(Day_ecdf$y), scale=diff(range(Day_ecdf$y)))
ED_ecdf$y <- scale(ED_ecdf$y, center=min(ED_ecdf$y), scale=diff(range(ED_ecdf$y)))
WD_ecdf$y <- scale(WD_ecdf$y, center=min(WD_ecdf$y), scale=diff(range(WD_ecdf$y)))


## take the inverse of a cdf, such that the p(x) > or = log wave event
NF_ecdf$inv_y <- ((NF_ecdf$y - max(NF_ecdf$y)) * (-1))
WEU_ecdf$inv_y <- ((WEU_ecdf$y - max(WEU_ecdf$y)) * (-1))
WEK_ecdf$inv_y <- ((WEK_ecdf$y - max(WEK_ecdf$y)) * (-1))
Day_ecdf$inv_y <- ((Day_ecdf$y - max(Day_ecdf$y)) * (-1))
ED_ecdf$inv_y <- ((ED_ecdf$y - max(ED_ecdf$y)) * (-1))
WD_ecdf$inv_y <- ((WD_ecdf$y - max(WD_ecdf$y)) * (-1))


## add site name to new data frames
NF_ecdf$Site <- "NavFac"
WEU_ecdf$Site <- "West End Urchin"
WEK_ecdf$Site <- "West End Kelp"
Day_ecdf$Site <- "Daytona"
ED_ecdf$Site <- "East Dutch"
WD_ecdf$Site <- "West Dutch"


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat_ecdf <- rbind(NF_ecdf, WEU_ecdf, WEK_ecdf, Day_ecdf, ED_ecdf, WD_ecdf)
dat_ecdf$Site <- factor(dat_ecdf$Site, levels=c("NavFac", "West End Urchin", "West End Kelp", 
                                                "Daytona", "East Dutch", "West Dutch"))


## plot all 
p1 <- ggplot(dat_ecdf, aes(x, inv_y, color=Site)) + my.theme +
  geom_line(lwd=1, alpha=.8) +
  scale_color_manual(values=pal_All) +
  xlab("total sheephead per site (per 0.1 hectare)") + ylab("inverse empirical CDF") +
  scale_x_continuous(labels=c("0","3.16","10","31.6","100","316.2")) 

print(p1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## temporal dynamics figure to superimpose above Sheephead time series ~~~~~~~~~
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat1 <- read.csv("NMDS_coordinates.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat1$SITE <- factor(dat1$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"), 
                   labels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                            "Daytona", "East Dutch", "West Dutch"))


## subset by site
NF1 <- filter(dat1, SITE %in% c("NavFac"))
WEK1 <- filter(dat1, SITE %in% c("WestEnd Kelp"))
Day1 <- filter(dat1, SITE %in% c("Daytona"))
ED1 <- filter(dat1, SITE %in% c("East Dutch"))
WEU1 <- filter(dat1, SITE %in% c("WestEnd Urchin"))
WD1 <- filter(dat1, SITE %in% c("West Dutch"))
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
graphics.off()
windows(h=4,w=10, record=TRUE)

## Plot Navfac 
p1 <- ggplot(NF1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_NF) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +  ylim(-1.25, 1.25) +
  xlim(2,85) +
  #ylab("NMDS Axis-1: System State") +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p1)





p2 <- ggplot(WEU1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WEU) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +  ylim(-1.25, 1.25) +
  ylab("NMDS Axis-1: System State") +
  xlim(2,85) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p2)




## Plot WestEnd Kelp 
p3 <- ggplot(WEK1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WEK) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +  ylim(-1.25, 1.25) +
  ylab("NMDS Axis-1: System State") +
  guides(color=FALSE) +
  xlim(2,85) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p3)




p4 <- ggplot(Day1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_Day) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +  ylim(-1.25, 1.25) +
  xlim(2,85) +
  ylab("NMDS Axis-1: System State") +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p4)





p5 <- ggplot(ED1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_ED) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +
  ylim(-1.25, 1.25) +
  xlim(2,85) +
  ylab("NMDS Axis-1: System State") +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p5)




p6 <- ggplot(WD1, aes(PERIOD,NMDS1,color=TRANSECT)) +
  #geom_vline(xintercept=c(0.5,-.75),size=.1, color="gray50") +
  scale_color_manual(values=pal_WD) +
  #geom_point() +
  geom_path() +
  theme_bw() +
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_rect(fill = "transparent", colour = NA),
        plot.background = element_rect(fill = "transparent", colour = NA)) +  ylim(-1.25, 1.25) +
  ylab("NMDS Axis-1: System State") +
  xlim(2,85) +
  guides(color=FALSE) +
  theme(axis.title.y = element_blank(),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank()) +
  annotate("text", x=80,y=.7,label="Barren",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-.15,label="Mixed",size=4.75, hjust=.2) +
  annotate("text", x=80,y=-1,label="Algae",size=4.75, hjust=.2) 
#annotate("text", x=0,y=72.5,label="NavFac",size=4.75,vjust=2.2,hjust=.5)
print(p6)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## CA sheephead time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat2 <- read.csv("CAsheephead.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat2$log <- log10(dat2$SumAdultCount+1)

dat2$Site <- factor(dat2$Station, 
                   levels=c("1","2","3","4","5","6"),
                   labels=c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))

dat2 <- filter(dat2, Site %in% c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))


NF2 <- filter(dat2, Site %in% c("NavFac"))
WEK2 <- filter(dat2, Site %in% c("West End Kelp"))
WEU2 <- filter(dat2, Site %in% c("West End Urchin"))
ED2 <- filter(dat2, Site %in% c("East Dutch"))
WD2 <- filter(dat2, Site %in% c("West Dutch"))
Day2 <- filter(dat2, Site %in% c("Daytona"))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=6,w=10, record = T)

splots.theme <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank()) 
  

## individual sites: CA sheephead abundance (natural = SumAdultCount or log scale)
s1 <- ggplot(NF2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "NavFac") +  scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme + theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines")) +  
  scale_x_continuous(limits=c(3,85)) 
print(s1)  


s2 <- ggplot(WEU2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "West End Urchin") +   scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme + scale_x_continuous(limits=c(3,85)) 
print(s2)


s3 <- ggplot(WEK2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "West End Kelp") +   scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme + scale_x_continuous(limits=c(3,85)) 
print(s3)


s4 <- ggplot(Day2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "Daytona") + scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme + theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines")) +  
  scale_x_continuous(limits=c(3,85)) 
print(s4)


s5 <- ggplot(ED2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "East Dutch") + scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme + theme(plot.margin = unit(c(.1,.1,2.5,.1), "lines")) +  
  scale_x_continuous(limits=c(3,85)) 
print(s5)


s6 <- ggplot(WD2, aes(x = Period, y = log)) + geom_line(lwd=1) +
  labs(title = "West Dutch") + scale_y_continuous(limits = c(0, 3)) + 
  my.theme + splots.theme +  theme(plot.margin = unit(c(.1,.1,2.5,.1), "lines")) +  
  scale_x_continuous(limits=c(3,85)) 
print(s6)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## combine time series with temporal dynamics ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
C1 <- s1 + annotation_custom(ggplotGrob(p1), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C1)

C2 <- s2 + annotation_custom(ggplotGrob(p2), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C2)

C3 <- s3 + annotation_custom(ggplotGrob(p3), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C3)

C4 <- s4 + annotation_custom(ggplotGrob(p4), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C4)

C5 <- s5 + annotation_custom(ggplotGrob(p5), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C5)

C6 <- s6 + annotation_custom(ggplotGrob(p6), xmin= -2, xmax=90, ymin= 1.8, ymax= 3.5)
print(C6)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## aggregate all plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
windows(h=6,w=16, record=TRUE)


p10 <- egg::ggarrange(tag_facet(C1 + facet_wrap(~"Period"), tag_pool="a"),
                 tag_facet(C2 + facet_wrap(~"Period"), tag_pool="b"),
                 tag_facet(C3 + facet_wrap(~"Period"), tag_pool="c"),
                 tag_facet(C4 + facet_wrap(~"Period"), tag_pool="d"),
                 tag_facet(C5 + facet_wrap(~"Period"), tag_pool="e"),
                 tag_facet(C6 + facet_wrap(~"Period"), tag_pool="f"),
                 nrow=2, ncol=3)

print(p10)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom axes -- started here and finished in adobe PDF ~~~~~~~~~~~~~~~~~~~~~~~
seq <- seq(1980, 2020, by = 4)

lab<-grid.xaxis(at=c(seq), 
                vp=vpStack(viewport(width=unit(20.4,"lines")),
                           viewport(y=1, xscale = c(1980,2018), just="left")))

lab2<-grid.yaxis(at=c(1,2,3,4), 
                vp=vpStack(viewport(width=unit(1,"lines")),
                           viewport(x=1, yscale = c(0,4), just="left")))

print(p10)

vp1 <- viewport(x = 0.065, y = 0.08, width = .9, height = .575, just="center")
vp2 <- viewport(x = .0475, y = 0.23, width = 1, height = .46, just="center")

pushViewport(vp1)
grid.draw(lab)
## save pdf and copy xaxis over to center and right column 

print(p10)
name<-grid.text(label="total CA sheephead per 0.1 hectare", x=.01,y=.5,rot=90)
pushViewport(vp2)
grid.draw(lab2)
## save pdf and copy yaxis to upper row
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
