setwd()

library(vegan)
library(betapart)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(adespatial)
library(fitdistrplus)

# NMDS ------------------------------------------------

#Read community matrix
abund<-read.csv2("commat.csv", sep=",")

set.seed(321)
abund.mds <- metaMDS(abund[,3:33],try = 999) ##creates a Bray Curtis dissimilarity/distance matrix
plot(abund.mds,type= "t")
plot(abund.mds,"sites")

#Using the scores function from vegan to extract the site scores and convert to a data.frame
abund.scores <- scores(abund.mds) %>% 
  as.data.frame() %>% 
  bind_cols(abund[,2])

colnames(abund.scores) <- c("NMDS1", "NMDS2", "Habitat")


my_nmds_plot <- ggplot(data = abund.scores,
                       aes(x = NMDS1,
                           y = NMDS2,
                           text = Habitat)) +
  #geom_text(aes(label = Habitat), size= 2, hjust = 0, nudge_x  = 0.04, check_overlap = TRUE) +
  geom_point(aes(color = Habitat),
             size = 4) + theme_classic() +
  scale_color_manual(labels = c("Open savanna", "Wooded savanna"), values = c("#D55E00", "#009E73")) + 
  theme(legend.title = element_text(size=12)) +
  theme(legend.text = element_text(size=12))
my_nmds_plot

# Nestedness and Turnover -------------------------------------------------

#betapart
dist<-bray.part(abund[,3:33])
bd <- beta.multi.abund(dist[[3]], index.family="bray")

# IVI ---------------------------------------------------------------------
library(BiodiversityR)
library(tidyverse)
library(dplyr)
rfivi <- read.csv("ivirf.csv")
clivi <- read.csv("ivicl.csv")

a <- rfivi %>% as.data.frame() %>% 
  mutate(count=rep(1, each= nrow(.))) %>%
importancevalue( site= "Sitecode", species ="Species", 
                  count= "count", basal= "basal", factor="", level= "") %>% head(10)

b<- clivi %>% as.data.frame() %>% 
  mutate(count=rep(1, each= nrow(.))) %>%
  importancevalue( site= "Sitecode", species ="Species", 
                   count= "count", basal= "basal", factor="", level= "") %>% head(10)


# Size structure ----------------------------------------------------------

sizes_g <- read.csv("sizes_open savanna.csv")
sizes_f <- read.csv("sizes_wooded savanna.csv")

Sizes_g <- sizes_g$Size
hist(Sizes_g, main="Open savanna",
     xlab="Size class in cm",
     xlim=c(0, 180),
     col="#D55E00",  freq=FALSE, breaks=25)
lines(density(Sizes_g, adjust=2, cut= 0), col = "blue", lwd = 2)


Sizes_f <- sizes_f$Size
hist(Sizes_f, main="Wooded savanna",
     xlab="Size class in cm",
     xlim=c(0, 180),
     col="#009E73",  freq=FALSE, breaks=25)
lines(density(Sizes_f, adjust = 2, cut = 0), col = "blue", lwd = 2)

fit_fw  <- fitdist(Sizes_f, "weibull")
summary(fit_fw)

fit_gw  <- fitdist(Sizes_g, "weibull")
summary(fit_gw)

par(mfrow=c(2,1))
plotdist(Sizes_f,demp = TRUE)

denscomp(list(fit_fw), main="Wooded savanna", breaks=15, xlab="Size class in cm")
denscomp(list(fit_gw), main="Open savanna", breaks=15, xlab="Size class in cm")
