## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Fig S1 a, b for SubstrateComplexity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(tidyverse)
library(dplyr)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
dat <- read.csv("NMDS_coordinates.csv", header = TRUE)

my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title = element_text(size=16),
                 axis.text = element_text(size=14),
                 plot.title = element_text(size=16),
                 legend.position="none") 
## END initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## plotting prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$rugosity <- dat$mean_rug / 10


## reorder transects: NF, WEK, WEU, Day, ED, WD
dat$ID <- factor(dat$ID, levels=c("1_10R", "1_22R", "1_32L", "1_39R", "1_45R",
                                  "3_10R", "3_22R", "3_32L", "3_39R", "3_45L",
                                  "2_10L", "2_22L", "2_32R", "2_39L", "2_45L",
                                  "6_10R", "6_22L", "6_22R", "6_32L", "6_39L",
                                  "5_10R", "5_22R", "5_32L", "5_39R", "5_45R",
                                  "4_10R", "4_22L", "4_32L", "4_39L", "4_45L"))


## custom color palette, ordered: NF, WEK, WEU, Day, ED, WD
cols <- scale_color_manual(values=c("#390b0b", "#721716", "#ab2221", "#cb5151", "#df9392", 
                                    "#301800", "#773b00", "#be5e00", "#f0841a", "#f5ad66", 
                                    "#141100", "#3d3400", "#7b6800", "#b99c00", "#e6d680", 
                                    "#131a0d", "#304221", "#4d6934", "#708f54", "#a0b58d", 
                                    "#003446", "#00536f", "#3386a2", "#66a4b9", "#99c3d1", 
                                    "#2b112b", "#552255", "#803280", "#a560a5", "#c79cc7")) 
## END plotting prep ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## custom functions and data wrangling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## calculate loess smoother 
loess.f <- function(y, x){
  t1 <- loess(y ~ x, span = 1)
  t2 <- predict(t1)
  out <- return(t2)
}


## create a data frame with single maximum velocity to merge later 
high.vel <- dat %>%
  group_by(ID) %>%
  filter(Vel==max(Vel)) %>%
  arrange(ID)


## select highest 5 velocities of community shift 
top.5 <- dat %>%
  dplyr::group_by(ID) %>%
  dplyr::arrange(Vel, .by_group=TRUE) %>%
  dplyr::top_n(5, Vel) %>%
  dplyr::summarize(avg = mean(Vel))


## join new data with previous table to provide metadata
joined <- merge(high.vel, top.5, by.y="ID")


## calculate loess-smoother 
loess_top.5 <- loess.f(joined$avg, joined$rugosity)
## END data wrangling ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## create and visualize plots ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## plot transects ordered by rugosity with velocity loess-smoothed 
p1 <- ggplot(joined, aes(x=rugosity, y=loess_top.5)) + my.theme +
  geom_line(color = "black", lwd = 1.5) + geom_line(color = "white", lwd = .75) +
  geom_point(data=joined, aes(x=rugosity, y=avg, color=ID)) + cols + 
  xlab("Rugosity") + ylab("Average of five highest velocities") 


## plot rugosity vs NMDS Axis-1: System State
p2 <- ggplot(dat, aes(x=rugosity, y=NMDS1, color=ID)) + my.theme + cols +
  geom_point(position=position_jitter(w=0.03), alpha=.5) +
  xlab("Rugosity") + ylab("NMDS Axis-1 System State") 
  

## visualize 
graphics.off()
windows(8,5,record=T)
print(p1)
print(p2)
## END graphics visualization ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
