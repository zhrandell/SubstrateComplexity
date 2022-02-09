## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Fig S6 from Substrate Complexity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

#library(MASS)
#library(fitdistrplus)
library(gridExtra)
#library(RColorBrewer)
library(ggpubr)
library(egg)
library(gtable)
library(grid)
library(magick)
library(magrittr)
library(here)
library(dplyr)
library(tidyverse)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat <- read.csv("CAsheephead.csv", header = TRUE)

my.theme <- theme(panel.grid.major = element_blank(), 
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
                   labels=c("NavFac", "West End Urchin","West End Kelp",
                            "West Dutch","East Dutch","Daytona"))


## once more to reorder the site labels (for plotting)
dat$Site <- factor(dat$Site, levels=c("NavFac","West End Urchin","West End Kelp","Daytona","East Dutch","West Dutch"))


## remove Sandy Cove obs
dat <- filter(dat, Site %in% c("NavFac","West End Urchin","West End Kelp",
                               "West Dutch","East Dutch","Daytona"))


## log transformation
dat$log <- log10(dat$SumAdultCount+1)


## filter by site
filter.site <- function(x){filter(dat, Site %in% c(x))}
NF <- filter.site("NavFac")
WEK <- filter.site("West End Kelp")
WEU <- filter.site("West End Urchin")
Day <- filter.site("Daytona")
ED <- filter.site("East Dutch")
WD <- filter.site("West Dutch")
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot kernal densities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom graphical params 
fill.cols <- scale_fill_manual(values=c("#ab2221", "#be5e00", "#7b6800", "#4d6934", "#3386a2", "#803280"))
cols <- scale_color_manual(values=c("#ab2221", "#be5e00", "#7b6800", "#4d6934", "#3386a2", "#803280"))
x.scale <- scale_x_continuous(labels=c("0","3.16","10","31.6","100","316.2")) 
no.legend.title <- theme(legend.title = element_blank())


## create plot 
p3 <- ggplot(dat, aes(log, fill=Site)) +
  geom_density(position="identity", color="black", alpha=0.3) +  
  my.theme + fill.cols + x.scale + no.legend.title + xlab("total sheephead per site (per 0.1 hectare)") 


## visualize 
graphics.off()
windows(w=8,h=5,record=TRUE)
print(p3)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## eCDFs~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ~~~~~~~~~~~~~~~~~~~~~~~~~~
## calcuate empirical cumulative density function and extract/sort/process values 
ecdf.f <- function(dat, site.name){
  t1 <- data.frame(x=unique(dat$log), y=ecdf(dat$log)(unique(dat$log))*length(dat$log))
  t2 <- scale(t1$y, center=min(t1$y), scale=diff(range(t1$y)))
  t3 <- ((t2 - max(t2)) * (-1))
  t4 <- cbind(t1, t2, t3)
  names(t4)[1:4] <- c("x", "unscaled_y", "y", "inv_y")
  t4$Site <- site.name
  out <- return(t4)
}


## apply ecdf.f() function to data
NF.e <- ecdf.f(NF, "NavFac")
WEK.e <- ecdf.f(WEK, "West End Kelp")
WEU.e <- ecdf.f(WEU, "West End Urchin")
Day.e <- ecdf.f(Day, "Daytona")
ED.e <- ecdf.f(ED, "East Dutch")
WD.e <- ecdf.f(WD, "West Dutch")


## combine cfd data frames into single df & and reorder sites in order to plot properly 
dat.e <- rbind(NF.e, WEK.e, WEU.e, Day.e, ED.e, WD.e)
dat.e$Site <- factor(dat.e$Site, levels=c("NavFac","West End Kelp","West End Urchin",
                                          "Daytona","East Dutch","West Dutch"))


## create plot and visualize
p1 <- ggplot(dat.e, aes(x, inv_y, color=Site)) +  
  geom_line(lwd=1, alpha=.8) + my.theme + cols + x.scale +
  xlab("total sheephead per site (per 0.1 hectare)") + ylab("inverse empirical CDF") 
  
print(p1)
## END eCDFs ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## temporal dynamics figure to superimpose above Sheephead time series ~~~~~~~~~
setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat.ts <- read.csv("NMDS_coordinates.csv", header = TRUE)


## site as.factor and ordered as desired for plotting 
dat.ts$SITE <- factor(dat.ts$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",
                                            "Daytona", "East Dutch", "West Dutch"))


## subset by site
filter.site.2 <- function(x){filter(dat.ts, SITE %in% c(x))}
NF.ts <- filter.site.2("NavFac")
WEK.ts <- filter.site.2("WestEnd Kelp")
WEU.ts <- filter.site.2("WestEnd Urchin")
Day.ts <- filter.site.2("Daytona")
ED.ts <- filter.site.2("East Dutch")
WD.ts <- filter.site.2("West Dutch")
## END data wrangling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plotting configuration for NMDS Axis-1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## custom colors 
NF.col <- scale_color_manual(values=c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392")) 
WEK.col <- scale_color_manual(values=c("#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66")) 
WEU.col <- scale_color_manual(values=c("#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680")) 
Day.col <- scale_color_manual(values=c("#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d")) 
ED.col <- scale_color_manual(values=c("#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1")) 
WD.col <- scale_color_manual(values=c("#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7")) 


## custom graphing commands
background <- theme(panel.border = element_blank(), 
                    panel.grid.major = element_blank(),
                    panel.grid.minor = element_blank(),
                    panel.background = element_rect(fill = "transparent", colour = NA),
                    plot.background = element_rect(fill = "transparent", colour = NA), 
                    legend.position = "none") 

no.axes <- theme(axis.title.y = element_blank(),
                 axis.title.x = element_blank(),
                 axis.text.x = element_blank(),
                 axis.text.y = element_blank(),
                 axis.ticks.x = element_blank(),
                 axis.ticks.y = element_blank()) 
  
x.lim <- xlim(2, 85)
y.lim <- ylim(-1.32, 1.2)


## add annotation designating system state
text.size <- 4.75
hjust <- 0.2
x1 <- 80
y1 <- 0.7
y2 <- -0.15
y3 <- -1.0

add.annotation <- function(label, y.state){
  annotate("text", x=x1, y=y.state, label=label, size=text.size, hjust=hjust)
  }

barren <- add.annotation("Barren", y1)
mixed <- add.annotation("Mixed", y2)
algae <- add.annotation("Algae", y3)
## END custom plotting params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Create plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.ts <- function(dat, cols){
  ggplot(dat, aes(PERIOD, NMDS1, color=TRANSECT)) + 
    cols + geom_path() + theme_bw() + background + x.lim + y.lim + no.axes +
    barren + mixed + algae
}


## call function to create plots
NF.p1 <- plot.ts(NF.ts, NF.col)
WEK.p1 <- plot.ts(WEK.ts, WEK.col)
WEU.p1 <- plot.ts(WEU.ts, WEU.col)
Day.p1 <- plot.ts(Day.ts, Day.col)  
ED.p1 <- plot.ts(ED.ts, ED.col)
WD.p1 <- plot.ts(WD.ts, WD.col)


## visualize 
graphics.off()
windows(h=4,w=10, record=TRUE)
print(NF.p1)
## END NMDS Axis-1 creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## CA sheephead time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat.sh <- read.csv("CAsheephead.csv", header = TRUE)


## as.factor 
dat.sh$Site <- factor(dat.sh$Station, 
                    levels=c("1","2","3","4","5","6"),
                    labels=c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))


## re-order for plotting 
dat.sh <- filter(dat.sh, Site %in% c("NavFac","West End Urchin","West End Kelp","West Dutch","East Dutch","Daytona"))


## log transform 
dat.sh$log <- log10(dat.sh$SumAdultCount+1)



filter.site.3 <- function(x){filter(dat.sh, Site %in% c(x))}
NF.sh <- filter.site.3("NavFac")
WEK.sh <- filter.site.3("West End Kelp")
WEU.sh <- filter.site.3("West End Urchin")
Day.sh <- filter.site.3("Daytona")
ED.sh <- filter.site.3("East Dutch")
WD.sh <- filter.site.3("West Dutch")
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plot sheephead time series ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
scale.x <- scale_x_continuous(limits=c(3, 85))
scale.y <- scale_y_continuous(limits=c(0, 3))


plot.sheephead <- function(dat, site.name){
  ggplot(dat, aes(x=Period, y=log)) + geom_line(lwd=1) + labs(title=site.name) + 
    scale.y + scale.x + no.axes + my.theme 
}

NF.p2 <- plot.sheephead(NF.sh, "NavFac") + theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines"))
WEK.p2 <- plot.sheephead(WEK.sh, "West End Kelp")
WEU.p2 <- plot.sheephead(WEU.sh, "West End Urchin")
Day.p2 <- plot.sheephead(Day.sh, "Daytona") + theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines")) 
ED.p2 <- plot.sheephead(ED.sh, "East Dutch") + theme(plot.margin = unit(c(.1,.1,2.5,.1), "lines"))
WD.p2 <- plot.sheephead(WD.sh, "West Dutch") + theme(plot.margin = unit(c(.1,.1,2.5,.1), "lines"))


## combine state dynamics along Axis-1 with sheephead time series
combine.plot <- function(f1, f2){
  f1 + annotation_custom(ggplotGrob(f2), xmin=-2, xmax=90, ymin=1.8, ymax=3.5)
}


NF.p3 <- combine.plot(NF.p2, NF.p1)
WEK.p3 <- combine.plot(WEK.p2, WEK.p1)
WEU.p3 <- combine.plot(WEU.p2, WEU.p1)
Day.p3 <- combine.plot(Day.p2, Day.p1)
ED.p3 <- combine.plot(ED.p2, ED.p1)
WD.p3 <- combine.plot(WD.p2, WD.p1)

print(NF.p3)


windows(h=6,w=16, record=TRUE)
figS7 <- ggarrange(NF.p3, WEK.p3, WEU.p3, Day.p3, ED.p3, WD.p3, nrow=2)
## END plot creation ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom axis creation . . . start here and finish in adobe ~~~~~~~~~~~~~~~~~~~
seq <- seq(1980, 2020, by = 4)

lab<-grid.xaxis(at=c(seq), 
                vp=vpStack(viewport(width=unit(20.4,"lines")),
                           viewport(y=1, xscale = c(1980,2018), just="left")))

lab2<-grid.yaxis(at=c(1,2,3,4), 
                 vp=vpStack(viewport(width=unit(1,"lines")),
                            viewport(x=1, yscale = c(0,4), just="left")))

print(figS7)

vp1 <- viewport(x = 0.065, y = 0.08, width = .9, height = .575, just="center")
vp2 <- viewport(x = .0475, y = 0.23, width = 1, height = .46, just="center")

pushViewport(vp1)
grid.draw(lab)
## save pdf and copy xaxis over to center and right column 

print(figS7)
name<-grid.text(label="total CA sheephead per 0.1 hectare", x=.01,y=.5,rot=90)
pushViewport(vp1)
pushViewport(vp2)
grid.draw(lab2)
## save pdf and copy yaxis to upper row ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
