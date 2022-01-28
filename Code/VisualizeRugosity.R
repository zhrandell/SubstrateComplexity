## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Visualize measurements of substrate rugosity ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## January 28th 2021 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

library(egg)
library(tidyverse)

setwd("D:/OneDrive/Active_Projects/SubstrateComplexity/Data")
relief_dat <- read.csv("SubstrateRugosity.csv", header = TRUE)
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## convert relief --> rugosity 
relief_dat$Rugosity <- relief_dat$RELIEF / 10


## calculate mean, sd, se 
new_rugosity <- ddply(relief_dat, c("SITE", "TRANSECT"), summarise, 
                      N = length(Rugosity),
                      mean_rug = mean(Rugosity),
                      sd_rug = sd(Rugosity),
                      se_rug = sd_rug / sqrt(N))


## filter data by site  
NF_rug <- filter(new_rugosity, SITE %in% c("NavFac"))
WEK_rug <- filter(new_rugosity, SITE %in% c("West End Kelp"))
Day_rug <- filter(new_rugosity, SITE %in% c("Daytona"))
ED_rug <- filter(new_rugosity, SITE %in% c("East Dutch"))
WEU_rug <- filter(new_rugosity, SITE %in% c("West End Urchin"))
WD_rug <- filter(new_rugosity, SITE %in% c("West Dutch"))


## reorder transects by ascending rugosity  
NF_rug$TRANSECT <- factor(NF_rug$TRANSECT, levels = NF_rug$TRANSECT[order(NF_rug$mean_rug)])
WEK_rug$TRANSECT <- factor(WEK_rug$TRANSECT, levels = WEK_rug$TRANSECT[order(WEK_rug$mean_rug)])
WEU_rug$TRANSECT <- factor(WEU_rug$TRANSECT, levels = WEU_rug$TRANSECT[order(WEU_rug$mean_rug)])
Day_rug$TRANSECT <- factor(Day_rug$TRANSECT, levels = Day_rug$TRANSECT[order(Day_rug$mean_rug)])
ED_rug$TRANSECT <- factor(ED_rug$TRANSECT, levels = ED_rug$TRANSECT[order(ED_rug$mean_rug)])
WD_rug$TRANSECT <- factor(WD_rug$TRANSECT, levels = WD_rug$TRANSECT[order(WD_rug$mean_rug)])
## END data configuration ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## graphical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pal.NF <- scale_fill_manual(values=c("#9f1009", "#5b0905", "#e3170d", "#e9453d", "#f18b86")) 
pal.WEK <- scale_fill_manual(values=c("#b34700", "#662900", "#ff6600", "#ff8533", "#ffb380"))
pal.WEU <- scale_fill_manual(values=c("#665c00", "#ffe600", "#b3a100", "#332e00", "#fff380"))
pal.Day <- scale_fill_manual(values=c("#91c591", "#0a2a0a", "#228b22", "#4ea24e", "#145314"))
pal.ED <- scale_fill_manual(values=c("#99c3d1", "#003446", "#3386a2", "#00536f", "#66a4b9")) 
pal.WD <- scale_fill_manual(values=c("#d6adeb", "#b870dc", "#4d1967", "#7a28a4", "#a347d2")) 

myTheme <- theme(panel.border = element_rect(color="gray50"), 
                 panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 axis.title.x = element_text(size=17),
                 axis.title.y = element_text(size=17),
                 axis.text.x = element_text(size=16), 
                 axis.text.y = element_text(size=16), 
                 plot.margin = margin(r=.5, l=.5, b=.5, t=.5, unit = "pt"), 
                 legend.position="none")

yScale <- scale_y_continuous(limits = c(.93,2.45), expand = c(0,0)) 

noYtitle <- theme(axis.title.y = element_blank(), 
                  axis.text.y = element_blank())

ptSize <- 3.5
ptShape <- 21
ptCol <- "black"
sdWidth <- 0.4
## END grahical params ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
graphics.off()
windows(h=4,w=24, record=TRUE)

## NavFac
p.NF <- ggplot(NF_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.NF + yScale + labs(x="NavFac") + labs(y="Rugosity") 
  
## West End Kelp 
p.WEK <- ggplot(WEK_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.WEK + yScale + noYtitle + labs(x="WestEnd Kelp") 

## West End Urchin 
p.WEU <- ggplot(WEU_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.WEU + yScale + noYtitle + labs(x="WestEnd Urchin") 

## Daytona
p.Day <- ggplot(Day_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.Day + yScale + noYtitle + labs(x="Daytona") 

## East Dutch 
p.ED <- ggplot(ED_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.ED + yScale + noYtitle + labs(x="East Dutch") 

## West Dutch 
p.WD <- ggplot(WD_rug, aes(SITE, mean_rug)) +
  geom_errorbar(aes(x=TRANSECT, ymin=mean_rug-se_rug, ymax=mean_rug+se_rug), width=sdWidth) +
  geom_point(aes(x=TRANSECT, y=mean_rug, fill=TRANSECT), size = ptSize, shape=ptShape, color=ptCol) +
  theme_bw() + myTheme + pal.WD + yScale + noYtitle + labs(x="West Dutch") 

## print
rugosity <- ggarrange(p.NF, p.WEK, p.WEU, p.Day, p.ED, p.WD, nrow=1)
## END plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~





## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## END Script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ 
