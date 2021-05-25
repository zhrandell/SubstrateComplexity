## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Visualize measurements of substrate complexity i.e. rugosity ~~~~~~~~~~~~~ ##
## Modified May 25th 2021; zhr ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## initiate ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
rm(list = ls())

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
## END Script ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##



