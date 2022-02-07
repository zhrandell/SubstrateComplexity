## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Script to produce Figure 3, Urchin Densities for Substrate Complexity ms ~ ##
## updated Feb 3 2022; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(MASS)
library(fitdistrplus)
library(gridExtra)
library(ggpubr)
library(egg)
library(gtable)
library(grid)
library(magick)
library(magrittr)
library(here)
library(pryr)
library(ggbeeswarm)
library(tidyverse)

dataLocation <- "D:/OneDrive/Active_Projects/SubstrateComplexity/Data"
figLocation <- "D:/OneDrive/Active_Projects/SubstrateComplexity/Figures/ms_figs"

setwd(dataLocation)
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$SITE <- factor(dat$SITE, levels=c("NavFac", "WestEnd Kelp", "WestEnd Urchin",  
                                      "Daytona", "East Dutch", "West Dutch"))


# calculate total urchin densities, algae densities (unused here), and ratios
dat$totalUrchin <- dat$StrPur + dat$MesFra


## function to subset by site
filter.site <- function(x, y){
  filter(x, SITE %in% c(y))
}


NF <- filter.site(dat, "NavFac")
WEK <- filter.site(dat, "WestEnd Kelp")
WEU <- filter.site(dat, "WestEnd Urchin")
Day <- filter.site(dat, "Daytona")
ED <- filter.site(dat, "East Dutch")
WD <- filter.site(dat, "West Dutch")
## end data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom graphing params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## set window size
no.border <- theme(panel.border = element_blank())
add.border <- theme(panel.border = element_rect(color="gray50")) 
x.scale <- scale_x_continuous(limits=c(0, 1))
y.scale <- scale_y_continuous(limits=c(0, 2.5), expand=c(0, 0))
margin <- theme(plot.margin = unit(c(.1,.1,.1,.1), "lines")) 

my.theme <- theme(panel.grid.major = element_blank(), 
                  panel.grid.minor = element_blank(), 
                  panel.background = element_rect(fill = "transparent",colour = NA), 
                  plot.background = element_rect(fill = "transparent",colour = NA),
                  plot.title = element_blank(),
                  axis.title.y = element_blank(),
                  axis.title.x = element_blank(),
                  axis.text.x = element_blank(),
                  axis.text.y = element_blank(),
                  axis.ticks = element_blank(),
                  legend.position = "none") 


## coordinates to match labels up with kernal density placement
plot.xmin <- 0; 
plot.xmax <- 1
algae.basin <- 0.15; 
mixed.basin <- 0.5; 
urchin.basin <- 0.85;
scaler <- 5.75; 
width <- (plot.xmax - plot.xmin)/scaler
label.y <- 0.045; 
title.x <- 0.4; 
title.y <- 0.95;
title.size <- 14;
text.size <- 14;
title.col <- "black";
hjust <- 0;


## create a site label (title)
site.name <- function(title, title.x){
  grobTree(text_grob(title, x=title.x, y=title.y, hjust=hjust, size = title.size, color = title.col))
}


## create a letter label for each state 
state.label <- function(label, x.position){
  grobTree(text_grob(label, x=x.position, y=label.y, hjust=hjust, size=text.size, color=title.col))
}


## call function to create site titles 
NF.title <- site.name("NavFac", title.x)
WEK.title <- site.name("West End Kelp", title.x-.15)
WEU.title <- site.name("West End Urchin", title.x-.175)
Day.title <- site.name("Daytona", title.x-.05)
ED.title <- site.name("East Dutch", title.x-.05)
WD.title <- site.name("West Dutch", title.x-.05)


## create state labels and adjust position 
A <- state.label("A", algae.basin + .02)
M <- state.label("M", mixed.basin - 0.02)  
B <- state.label("B", urchin.basin - .04)


## x.axis position for black/red mean bars 
x.1 <- 0.09
x.2 <- 0.438
x.3 <- 0.79
## END custom graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## extract "peaks and troughs" from a kernal density plot with two states
extract.2 <- function(dat, key){
  x_dat <- density(dat)$x
  y_dat <- density(dat)$y
  peak_1_y <- which.max(y_dat)
  peak_1_x <- x_dat[peak_1_y]
  peak_2 <- max(y_dat[x_dat < key])
  peak_2_y <- which(y_dat == peak_2)
  peak_2_x <- x_dat[peak_2_y]
  valley_1 <- min(y_dat[x_dat < peak_1_x & x_dat > peak_2_x])
  valley_1_y <- which(y_dat == valley_1)
  valley_1_x <- x_dat[valley_1_y]
  out <- as.numeric(list(peak_1_x, peak_2_x, valley_1_x))
  return(out)
}


## this extract.2.V2 is a slightly modified version necessary for Daytona's altered kernal density structure 
extract.2.V2 <- function(dat, key){
  x_dat <- density(dat)$x
  y_dat <- density(dat)$y
  peak_1_y <- which.max(y_dat)
  peak_1_x <- x_dat[peak_1_y]
  peak_2 <- max(y_dat[x_dat > key]) ## switch to x_dat > key from x_dat < key
  peak_2_y <- which(y_dat == peak_2)
  peak_2_x <- x_dat[peak_2_y]
  valley_1 <- min(y_dat[x_dat > peak_1_x & x_dat < peak_2_x]) ## switch direction of both <, >
  valley_1_y <- which(y_dat == valley_1)
  valley_1_x <- x_dat[valley_1_y]
  out <- as.numeric(list(peak_1_x, peak_2_x, valley_1_x))
  return(out)
}


## extract "peaks and troughs" from a kernal density plot with three states 
extract.3 <- function(site, column, key1, key2){
  kd <- ggplot(site, aes(column)) +
    geom_density(adjust=0.8)
  kd.info <- ggplot_build(kd)
  df <- as.data.frame(kd.info[[1]])
  x_dat <- df$x
  y_dat <- df$y
  peak_y1 <- which.max(y_dat)
  peak_x1 <- x_dat[peak_y1]
  peak_2 <- max(y_dat[x_dat < key1])
  peak_y2 <- which(y_dat == peak_2)
  peak_x2 <- x_dat[peak_y2]
  valley_1 <- min(y_dat[x_dat < peak_x1 & x_dat > peak_x2])
  valley_y1 <- which(y_dat == valley_1)
  valley_x1 <- x_dat[valley_y1]
  peak_3 <- max(y_dat[x_dat > key2]) 
  peak_y3 <- which(y_dat == peak_3)
  peak_x3 <- x_dat[peak_y3]
  valley_2 <- min(y_dat[x_dat > peak_x1 & x_dat < peak_x3])
  valley_y2 <- which(y_dat == valley_2)
  valley_x2 <- x_dat[valley_y2]
  out <- as.numeric(list(peak_x1, peak_x2, valley_x1, peak_x3, valley_x2)) 
  return(out)
}


## check and make sure kernal density max / mins were extracted properly (used for diagnostic purposes)
plot.borders <- function(df, dat, list){
  ggplot(df, aes(dat)) + geom_density(data=df, aes(x=dat)) +
    geom_vline(xintercept = borders[1]) + 
    geom_vline(xintercept = borders[2]) +
    geom_vline(xintercept = borders[3])
}


## function to plot geom_beeswarm of total urchins -- creates core unit of final plot 
y.lim <- ylim(0, 2632)
proxy.col <- "white"
proxy.alph <- 0
pt.type <- 21
pt.fill <- "#C2C2C2"
pt.col <- "black"
pt.size <- 0.8
cex <- 5
stroke <- 0.35
ptStyle <- "random"

beeswarm.f <- function(dat){
  ggplot(dat) +
    geom_boxplot(data=dat, aes(y=totalUrchin), fill=proxy.col, alpha=proxy.alph, color=proxy.col, outlier.shape=NA) + 
    geom_beeswarm(data=dat, aes(x=0, y=totalUrchin), shape=pt.type, fill=pt.fill, color=pt.col, size=pt.size, groupOnX=NULL, cex=cex, stroke=stroke, priority=c(ptStyle)) +
    theme_bw() + no.border + my.theme + y.lim #+ scale_y_continuous(breaks=c(0, 650, 1300, 1950, 2600))
}


## extract median urchin density from geom_beeswarm() for two states
urchin.mean.2 <- function(plot1, plot2){
  mean1 <- ggplot_build(plot1)$data[[1]][[3]]
  mean2 <- ggplot_build(plot2)$data[[1]][[3]]
  out <- as.numeric(list(mean1, mean2))
}


## extract median urchin density from geom_beeswarm() for three states
urchin.mean.3 <- function(plot1, plot2, plot3){
  mean1 <- ggplot_build(plot1)$data[[1]][[3]]
  mean2 <- ggplot_build(plot2)$data[[1]][[3]]
  mean3 <- ggplot_build(plot3)$data[[1]][[3]]
  out <- as.numeric(list(mean1, mean2, mean3))
}


## add single pt onto beeswarm ID'ing location of median 
mean.pts <- function(plot, mu){
  plot <- plot + geom_point(aes(x=0, y=mu), color="cyan", size=0.5)
}


## add black bar as an outline around final mean plotted on aggregated beeswarm figure
bar.width <- 0.14
outline <- 0.0025

black.bar <- function(x, y){
  annotate("segment", x=(x - outline), xend=(x + bar.width + outline), y=y, yend=y, size=1.5, color="black") 
}


## add red bar as the primary geom_segment highlighting median urchin density 
red.bar <- function(x, y){
  annotate("segment", x=(x), xend=(x + bar.width), y=y, yend=y, size=1, color="red") 
}


## function that invokes all previous functions to create beewarm plot; 2 states
apply.2 <- function(dat, site, column, state1, state2, key){
  borders <- extract.2(dat, key)
  state1 <- filter(site, column > borders[3])
  state2 <- filter(site, column < borders[3])
  p1 <- beeswarm.f(state1)
  p2 <- beeswarm.f(state2)
  means <- urchin.mean.2(p1, p2)
  p1 <- mean.pts(p1, means[1])
  p2 <- mean.pts(p2, means[2])
  return(list(p1, p2))
}


## function that invokes all previous functions to create beewwarm plot (Daytona alternative); 2 states
apply.2.V2 <- function(dat, site, column, state1, state2, key){
  borders <- extract.2.V2(dat, key)
  state1 <- filter(site, column > borders[3])
  state2 <- filter(site, column < borders[3])
  p1 <- beeswarm.f(state1)
  p2 <- beeswarm.f(state2)
  means <- urchin.mean.2(p1, p2)
  p1 <- mean.pts(p1, means[1])
  p2 <- mean.pts(p2, means[2])
  return(list(p1, p2))
}


## function that invokes all previous functions to create beeswarm plot; 3 states
apply.3 <- function(site, column, key1, key2, urchin.state, mixed.state, algae.state){
  borders <- extract.3(site, column, key1, key2)
  urchin.state <- filter(site, column > borders[5])
  mixed.state <- filter(site, column < borders[5] & column > borders[3])
  algae.state <- filter(site, column < borders[3])
  p1 <- beeswarm.f(urchin.state)
  p2 <- beeswarm.f(mixed.state)
  p3 <- beeswarm.f(algae.state)
  means <- urchin.mean.3(p1, p2, p3)
  p1 <- mean.pts(p1, means[1])
  p2 <- mean.pts(p2, means[2])
  p3 <- mean.pts(p3, means[3])
  return(list(p1, p2, p3))
}
## END functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## invoke functions to extract borders, calculate data frame, create figs ~~~~~~
NF_figs <- apply.2(NF$NMDS1, NF, NF$NMDS1, NF_urchin, NF_mixed, 0)
ED_figs <- apply.2(ED$NMDS1, ED, ED$NMDS1, ED_urchin, ED_mixed, -0.5)
Day_figs <- apply.2.V2(Day$NMDS1, Day, Day$NMDS1, Day_urchin, Day_mixed, 0.5)
WEK_figs <- apply.3(WEK, WEK$NMDS1, -0.75, 0.5, WEK_urchin, WEK_mixed, WEK_algae)
WEU_figs <- apply.3(WEU, WEU$NMDS1, -0.75, 0.5, WEU_urchin, WEU_mixed, WEU_algae)
## END calculations ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




lab.x <- .04
lab.y <- 2.35
lab.size <- 8


fig.panel <- function(x){
  annotate("text", x = lab.x, lab.y, label=x, size = lab.size)
}


lab.A <- fig.panel("A")
lab.B <- fig.panel("B")
lab.C <- fig.panel("C")
lab.D <- fig.panel("D")
lab.E <- fig.panel("E")
lab.F <- fig.panel("F")




## aggregate all information into final site-level figures ~~~~~~~~~~~~~~~~~~~~~
## 1 state (West Dutch only)
## calculate relevant information (extract.X() not required for single basin)
p1 <- beeswarm.f(WD)
means <- urchin.mean.2(p1, p1)
p1 <- mean.pts(p1, means[1])


## plot West Dutch 
y.2 <- 0.355

p.WD <- ggplot(WD, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + margin + 
  annotation_custom(ggplotGrob(p1), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(WD.title) + annotation_custom(M) + black.bar(x.2, y.2) + red.bar(x.2, y.2) + lab.F
#print(p.WD)


## 2 states: NavFac, East Dutch, Daytona 
## plot NavFac 
y.2 <- 0.325
y.3 <- 0.71

p.NF <- ggplot(NF, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + 
  theme(plot.margin = unit(c(.1,.1,.1,4.5), "lines")) +
  annotation_custom(ggplotGrob(NF_figs[1][[1]]), xmin = (urchin.basin - width), xmax = (urchin.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(NF_figs[2][[1]]), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(NF.title) + annotation_custom(M) + annotation_custom(B) +
  black.bar(x.2, y.2) + red.bar(x.2, y.2) + black.bar(x.3, y.3) + red.bar(x.3, y.3) + lab.A
print(p.NF)


## plot East Dutch 
y.1 <- 0.21
y.2 <- 0.35

p.ED <- ggplot(ED, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + margin +
  annotation_custom(ggplotGrob(ED_figs[1][[1]]), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(ED_figs[2][[1]]), xmin = (algae.basin - width), xmax = (algae.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ED.title) + annotation_custom(M) + annotation_custom(A) +
  black.bar(x.1, y.1) + red.bar(x.1, y.1) + black.bar(x.2, y.2) + red.bar(x.2, y.2) + lab.E
print(p.ED)


## plot Daytona 
y.2 <- 0.57
y.3 <- 0.825

p.Day <- ggplot(Day, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + margin + 
  annotation_custom(ggplotGrob(Day_figs[1][[1]]), xmin = (urchin.basin - width), xmax = (urchin.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(Day_figs[2][[1]]), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(Day.title) + annotation_custom(M) + annotation_custom(B) +
  black.bar(x.2, y.2) + red.bar(x.2, y.2) + black.bar(x.3, y.3) + red.bar(x.3, y.3) + lab.D
#print(p.Day)


## 3 states: West End Urchin, West End Kelp 
## plot West End Kelp 
y.1 <- 0.20575
y.2 <- 0.31
y.3 <- 0.685

p.WEK <- ggplot(WEK, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + margin + 
  annotation_custom(ggplotGrob(WEK_figs[1][[1]]), xmin = (urchin.basin - width), xmax = (urchin.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(WEK_figs[2][[1]]), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(WEK_figs[3][[1]]), xmin = (algae.basin - width), xmax = (algae.basin + width), ymin = 0, ymax = 2.5) +
  black.bar(x.1, y.1) + red.bar(x.1, y.1) + black.bar(x.2, y.2) + red.bar(x.2, y.2) + black.bar(x.3, y.3) + red.bar(x.3, y.3) + 
  annotation_custom(WEK.title) + annotation_custom(A) + annotation_custom(M) + annotation_custom(B) + lab.B
#print(p.WEK)


## plot West End Urchin
y.1 <- 0.20575
y.2 <- 0.365
y.3 <- 0.85

p.WEU <- ggplot(WEU, aes(NMDS1)) + 
  add.border + theme_bw() + x.scale + y.scale + my.theme + margin + 
  annotation_custom(ggplotGrob(WEU_figs[1][[1]]), xmin = (urchin.basin - width), xmax = (urchin.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(WEU_figs[2][[1]]), xmin = (mixed.basin - width), xmax = (mixed.basin + width), ymin = 0, ymax = 2.5) +
  annotation_custom(ggplotGrob(WEU_figs[3][[1]]), xmin = (algae.basin - width), xmax = (algae.basin + width), ymin = 0, ymax = 2.5) +
  black.bar(x.1, y.1) + red.bar(x.1, y.1) + black.bar(x.2, y.2) + red.bar(x.2, y.2) + black.bar(x.3, y.3) + red.bar(x.3, y.3) + 
  annotation_custom(WEU.title) + annotation_custom(A) + annotation_custom(M) + annotation_custom(B) + lab.C
#print(p.WEU)
## END final arrangement of individual site-level figures ~~~~~~~~~~~~~~~~~~~~~~





## plot all six figures ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=3,w=15, record=TRUE)


## custom axis for later call 
new<-grid.yaxis(at=c(0,650,1300,1950,2600), 
                vp=vpStack(viewport(width=unit(2,"lines")),
                           viewport(x=1, yscale = c(-130,1700), just="left")))


## aggregate all plot 
#p2 <- ggarrange(tag_facet(p.NF + facet_wrap(~"NMDS1"), tag_pool = "a"),
#                tag_facet(p.WEK + facet_wrap(~"NMDS1"), tag_pool = "b"),
#                tag_facet(p.WEU + facet_wrap(~"NMDS1"), tag_pool = "c"),
#                tag_facet(p.Day + facet_wrap(~"NMDS1"), tag_pool = "d"),
#                tag_facet(p.ED + facet_wrap(~"NMDS1"), tag_pool = "e" ),
#                tag_facet(p.WD + facet_wrap(~"NMDS1"), tag_pool = "f" ),
#                nrow=1)



fig4 <- ggarrange(p.NF, p.WEK, p.WEU, p.Day, p.ED, p.WD, nrow=1)

## add axis legend 
name<-grid.text(label="total Urchin abundance", x=.01,y=.5,rot=90)

## add previously created axis
sample_vp <- viewport(x = .545, y = 0.3475, width = 1, height = .575, just = c("right", "center"))
pushViewport(sample_vp)
grid.draw(new)
grid.draw(name)
popViewport()




setwd(figLocation)




ggplot2::ggsave(filename="fig4-urchdens.pdf", 
       plot=fig4,
       device=cairo_pdf, 
       height=3,
       width=15,
       units="in", 
       useDingbats=FALSE)


## END of final figure~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

print(fig4)



## EXTRA code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## This code walks through the steps that's otherwise wrapped within the apply.f() function
## running this code allows one to double check the placement of "peaks" and "valleys", 
## the max and min pts along a kernal density, identified with extract.2() or extract.3(), for 
## 2 or 3 basins of attraction, respectively
## uncomment the lines with a single "#"

## apply functions to extract kernal density peaks/troughs and plot 
## apply function and store relevant output 
#borders <- extract.2(NF$NMDS1, 0)


## apply function and examine plot 
#graphics.off()
#windows(h=5,w=5, record=TRUE)


#a1 <- plot.borders(df = NF, dat = NF$NMDS1, list = borders)
#print(a1)


## create new data frames containing data for the barren and mixed states 
#NF_urchin <- filter(NF, NMDS1 > borders[3])
#NF_mixed <- filter(NF, NMDS1 < borders[3])


## create beeswarm cluster from the mixed state
#p1 <- beeswarm.f(NF_mixed)
#p2 <- beeswarm.f(NF_urchin)


## extract means
#means <- urchin.mean.2(p1, p2)


## add mean lines 
#p1 <- mean.pts(p1, means[1])
#p2 <- mean.pts(p2, means[2])


## check plots
#print(p1)
#print(p2)
## END extra code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END of script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
