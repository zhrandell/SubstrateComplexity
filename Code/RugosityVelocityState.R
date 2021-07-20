## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Velocity/Rugosity/State SOM figs for SubstrateComplexity ms ~~~~~~~~~~~~~~ ##
## Modified July 16th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
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
library(plyr)

## ordination data
setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

## rugosity data 
#setwd("D:/OneDrive/Active_Projects/Substrate_Complexity/Data")
#relief_dat <- read.csv("SubstrateRugosity.csv", header = TRUE)

## external plotting window 
graphics.off()
windows(w=8,h=5,record=TRUE)

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16)) 
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## process rugosity values and add to ordination data sheet 
dat$rugosity <- dat$mean_rug / 10


## reorder transects: NF, WEK, WEU, Day, ED, WD
dat$ID <- factor(dat$ID, levels=c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                        "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                        "2_10L", "2_22L", "2_32R", "2_39L", "2_45L",
                                        "6_10R", "6_22L", "6_22R", "6_32L", "6_39L",
                                        "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                        "4_10R", "4_22L", "4_32L", "4_39L", "4_45L"))


# Custom color palette, ordered: NF, WEK, WEU, Day, ED, WD
pal_All <- c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392", #BE2625 used for tint and shade creation
             "#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66", #FF6600 ""  ""
             "#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680", #CDAD00 ""  ""
             "#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d", #608341 ""  ""
             "#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1", #00688B ""  ""
             "#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7") #9932CD ""  ""
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Velocity vs Rugosity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## select Max Velocity event 
highvel <- dat %>%
  group_by(ID) %>%
  filter(Vel==max(Vel)) %>%
  arrange(ID)


## loess smoother through max velocities (1 per transect)
smooth <- loess(highvel$Vel ~ highvel$rugosity, data = highvel, span=1)
loess <- predict(smooth)


## plot loess and raw velocities
c1 <- ggplot(highvel, aes(x=rugosity, y=loess)) + my.theme +
  geom_line(color = "black", lwd = 1.5) + geom_line(color = "white", lwd = .75) +
  geom_point(data=highvel, aes(x=rugosity, y=Vel, color=ID)) +
  scale_color_manual(values=pal_All) +
  xlab("Rugosity") + ylab("Max velocity") +
  theme(legend.position = "none")

print(c1)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## select upper X number of velocities ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
topX <- dat %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(Vel, .by_group=TRUE) %>%
  dplyr::top_n(5, Vel) %>%
  dplyr::summarize(mn = mean(Vel))


## join new data with previous table
joined <- merge(highvel, topX, by.y="ID")


## add loess smoother
smooth_upper <- loess(joined$mn ~ joined$rugosity, data = joined, span=1)
loess_upper <- predict(smooth_upper)


## plot loess and raw velocities
c1 <- ggplot(joined, aes(x=rugosity, y=loess_upper)) + my.theme +
  geom_line(color = "black", lwd = 1.5) + geom_line(color = "white", lwd = .75) +
  geom_point(data=joined, aes(x=rugosity, y=mn, color=ID)) +
  scale_color_manual(values=pal_All) +
  xlab("Rugosity") + ylab("Average of five highest velocities") +
  theme(legend.position = "none")
print(c1)
## END Velocity vs Rugosity plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Rugosity vs system state ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# plot 
c2 <- ggplot(dat, aes(x=rugosity, y=NMDS1, color=ID)) + my.theme +
  geom_point(position=position_jitter(w=0.03), alpha=.5) +
  scale_color_manual(values=pal_All) +
  xlab("Rugosity") + ylab("NMDS Axis-1 System State") +
  theme(legend.position = "none")
print(c2)
## end script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

