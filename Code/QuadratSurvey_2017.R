## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## Visualization / Analyses of 2017 Quadrat Survey around San Nicolas Island 
## zhr -- May 23rd 2021 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ ##





## startup ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## clear Global Environment 
rm(list = ls())

## load packages (if packages not installed: install.packages("librarynamehere"))
library(ggplot2)
library(dplyr)
library(vegan)
library(MASS)
library(fitdistrplus)
library(pryr)
library(ggbeeswarm)
library(tidyverse)
library(ggpubr)
library(beeswarm)

## close out window tab & open custom window
graphics.off()
windows(h=6,w=9, record=TRUE)

## set working directory
setwd("D:/Active_Projects/Substrate_Complexity/Code")

## load data frame
dat <- read.csv("QuadSurvey2017.csv", header = TRUE)

## set up custom ggplot theme 
my.theme = theme(panel.grid.major = element_blank(), 
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(), 
                 axis.line = element_line(colour = "black"),
                 axis.title=element_text(size=16),
                 axis.text=element_text(size=14),
                 plot.title = element_text(size=16))
## ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~




dat <- filter(dat, Zone %in% c("1","2","3")) #,"5","6"))


## Create new column with sites broken down by relief ~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$Ind <- NA
dat$Ind[dat$Relief=="Low" & dat$Site=="NavFac"]<-"NavFac_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="NavFac"]<-"NavFac_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="East Dutch"]<-"EastDutch_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="East Dutch"]<-"EastDutch_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="West Dutch"]<-"WestDutch_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="West Dutch"]<-"WestDutch_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="West End"]<-"WestEnd_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="West End"]<-"WestEnd_High"
dat$Ind[dat$Relief=="Low" & dat$Site=="Daytona"]<-"Daytona_Low"
dat$Ind[dat$Relief=="High" & dat$Site=="Daytona"]<-"Daytona_High"







## Configure data ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
dat$Ind <- factor(dat$Ind, levels = c("NavFac_Low", "NavFac_High", 
                                        "WestEnd_Low", "WestEnd_High",
                                        "Daytona_Low", "Daytona_High",
                                        "EastDutch_Low", "EastDutch_High",
                                        "WestDutch_Low", "WestDutch_High"))



#dat$Relief <- factor(dat$Relief, levels = c("Low","High"))
color_pal <- c("#9D1309","#0E808B")



## split up data types; transform coverage to 0-1 scale
metaData <- dat[,c(1:9,31)]
percent <- dat[,10:20] / 100
abundance <- dat[,21:30]
dat <- cbind(metaData,percent,abundance)




## close out window tab & open custom window
graphics.off()
windows(h=10,w=15, record=TRUE)



par(mfrow=c(2,2))





par(mfrow = c(2, 2), mar = c(2, 2, 0, 0))
for(i in 1:4) {
  beeswarm(Cucumber ~ Ind, data = dat, main = '(a)',
           cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
           method="hex", corral = "wrap",
           xlab="", ylab="Cucumber % cover", 
           labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
  legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
  bxplot(Cucumber ~ Ind, data = dat, probs=.5, add = TRUE)
  
  beeswarm(FleshyRed ~ Ind, data = dat, main = '(b)',
           cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
           method="hex", corral = "wrap",
           xlab="", ylab="FleshyRed % cover", 
           labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
  legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
  bxplot(FleshyRed ~ Ind, data = dat, probs=.5, add = TRUE)
  
  
  
  beeswarm(ArtCor ~ Ind, data = dat, main = '(c)',
           cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
           method="hex", corral = "wrap",
           xlab="", ylab="Articulated Coralline % cover", 
           labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
  legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
  bxplot(ArtCor ~ Ind, data = dat, probs=.5, add = TRUE)
  
  
  
  beeswarm(EncrustCor ~ Ind, data = dat, main = '(d)',
           cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
           method="hex", corral = "wrap",
           xlab="", ylab="Encrusting Coralline % cover", 
           labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
  legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
  bxplot(EncrustCor ~ Ind, data = dat, probs=.5, add = TRUE)
  
  
  
  mtext(paste0("(", letters[i], ")"), side = 3, adj = 0.05, 
        line = -1.3)
}




beeswarm(Cucumber ~ Ind, data = dat, main = '(a)',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Cucumber % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(Cucumber ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(FleshyRed ~ Ind, data = dat, main = '(b)',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="FleshyRed % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(FleshyRed ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(ArtCor ~ Ind, data = dat, main = '(c)',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Articulated Coralline % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(ArtCor ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(EncrustCor ~ Ind, data = dat, main = '(d)',
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Encrusting Coralline % cover", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(EncrustCor ~ Ind, data = dat, probs=.5, add = TRUE)







beeswarm(Cystoseira ~ Ind, data = dat,
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Cystoseira abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topright", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(Cystoseira ~ Ind, data = dat, probs=.5, add = TRUE)




beeswarm(LamSpp ~ Ind, data = dat,
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Laminaria spp. abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(LamSpp ~ Ind, data = dat, probs=.5, add = TRUE)




beeswarm(PurpUrch ~ Ind, data = dat,
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Purple Urchin abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(PurpUrch ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(RedUrch ~ Ind, data = dat,
         cex=1.3, pch=21, col = c(1,4), bg = "#00000020", 
         method="hex", corral = "wrap",
         xlab="", ylab="Red Urchin abundance", 
         labels=c("NavFac","","WestEnd","","Daytona","","EastDutch","","WestDutch"))
legend("topleft", legend=c("Low","High"),title="Rugosity",pch=21,col=c(1,4),cex=1.3)
bxplot(RedUrch ~ Ind, data = dat, probs=.5, add = TRUE)










## group by Site, Relief
bySite <- dat2 %>%
  group_by(Site, Relief) 

## set up empty columns for long form data
name <- c("sppName", "Count")
dat2[,name] <- NA

## convert to long form
longForm <- bySite %>%
  gather("sppName","Count",
  FleshyRed,EncrustCor,ArtCor,Dictyota,Sargassum,Sponge,Cucumber,Infaunal,Sessile,BareRock,Sand,
  PurpUrch,RedUrch,Megastraea,GiantKelp,Eisenia,Pteryogophora,Cystoseira,LamSpp,YoungLam,Desmarestia) %>%
  arrange(desc(sppName, count)) 


meanList <- longForm %>%
  group_by(Site, Relief, sppName) %>%
  summarize(avg = mean(Count))

meanList %>%
  group_by(Site, Relief)



## Plot ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
## Verticality
p1 <- ggplot(dat, aes(Site, Verticality, group = Relief, color = Relief)) +
  ylim(c(0,10)) +
  scale_color_manual(values=color_pal) +
  geom_beeswarm(dodge.width = .5, cex=.45, stroke=.35, priority = c("random"), groupOnX = TRUE) +
  my.theme
print(p1)




p2 <- ggplot(dat3, aes(Site, Cucumber, color = Relief, alpha=.5)) +
  ylim(c(0,1)) +
  scale_color_manual(values=color_pal) +
  geom_beeswarm(dodge.width = 1, cex=1, stroke=1, priority = c("random"), groupOnX = TRUE) +
  my.theme 
print(p2)




p2 <- ggplot(dat2, aes(Site, Cucumber, group = Relief, color = Relief, alpha=.5)) +
  ylim(c(0,1)) +
  scale_color_manual(values=color_pal) +
  geom_beeswarm(#dodge.width = .5, 
                cex=2, stroke=1, priority = c("random"), groupOnX = TRUE) +
  my.theme +
  stat_summary(fun = median, geom="crossbar", width=0.5)



print(p2)



p9 <- ggplot(dat3, aes(Site, Cucumber, color = Relief, alpha = .2)) + 
  my.theme +
  scale_color_manual(values=color_pal) +
  geom_boxplot() +
  geom_dotplot(fill="white", position = "dodge", 
               binaxis = "y",
               binwidth = 0.01,
               stackdir = "center") 

  

print(p9)









p9 <- ggplot(dat3, aes(Site, Cucumber, fill = Relief)) + 
  my.theme +
  scale_fill_manual(values=color_pal) +
  #geom_boxplot() +
  geom_dotplot(dotsize=.25,
               position=position_jitterdodge(jitter.height=0.03, jitter.width=0.01, dodge.width=1),  
               binaxis = "y",
               binwidth = .1,
               stackdir = "center",
               stackratio = 1, 
               method="dotdensity",
               binpositions = "all") 


print(p9)









p10 <- ggplot(dat3, aes(Site, Cucumber, fill = Relief)) + 
  my.theme +
  scale_fill_manual(values=color_pal) +
  geom_boxplot() 

print(p10)









range(dat2$Cucumber)




ggdotplot(dat3, x="Site", y="Cucumber", fill="Relief", 
          add="mean")




p3 <- ggplot(dat2, aes(Site, Cucumber, group = Relief, color=Relief)) +
  stat_summary(fun="mean", geom="point")
print(p3)



fig_dat1 <- ggplot_build(p2)$data[[1]]
fig_dat2 <- ggplot_build(p3)$data[[1]]






beeswarm(Cucumber ~ Ind, data = dat,
         cex=1.5, pch=21, col = 1:2, bg = "#00000050", 
         method="hex",
         corral = "wrap",
         side=0)
bxplot(Cucumber ~ Ind, data = dat, probs=.5, add = TRUE)



beeswarm(Cucumber ~ Site, data = dat2, pch=NA)

x0 <- dat$Cucumber[dat$Relief == "Low"]
x1 <- dat$Cucumber[dat$Relief == "High"]

beeswarm(x0, pch=16, method='center', side=1, col=1, at=1 + xinch(0.04), add = TRUE)
beeswarm(x1, pch=16, method='center', side=-1, col=2, at=1 - xinch(0.04), add = TRUE)


















# [1]
## Generate some random data
set.seed(123)
distro <- list(runif = runif(100, min = -3, max = 3), 
               rnorm = rnorm(100))
## Compare the five "corral" methods 
for (ii in c("none", "gutter", "wrap", "random", "omit"))
{
  beeswarm(distro, 
           pch = 21, col = 2:4, bg = "#00000050",
           corral = ii, 
           main = paste('corral = "', ii, '"', sep = ''))
}

# [2]
distributions <- list(runif = runif(100, min = -3, max = 3), 
                      rnorm = rnorm(100),
                      rlnorm = rlnorm(100, sdlog = 0.5))
## Demonstrate the 'corral' methods
par(mfrow = c(2,3))
beeswarm(distributions, col = 2:4, 
         main = 'corral = "none" (default)')
beeswarm(distributions, col = 2:4, corral = "gutter", 
         main = 'corral = "gutter"')
beeswarm(distributions, col = 2:4, corral = "wrap", 
         main = 'corral = "wrap"')
beeswarm(distributions, col = 2:4, corral = "random", 
         main = 'corral = "random"')
beeswarm(distributions, col = 2:4, corral = "omit", 
         main = 'corral = "omit"')
