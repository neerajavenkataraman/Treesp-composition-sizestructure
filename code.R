setwd("C:/Users/vn38fiva/Desktop/chapter 3/chap3analysis/chap3")

library(vegan)
library(betapart)
library(dplyr)
library(tidyr)
library(ggplot2)
library(forcats)
library(adespatial)
library(fitdistrplus)


# Make community matrix ---------------------------------------------------

#make a community matrix, plotwise, sitewise and landusewise
dt <- read.csv("plotwise.csv")
dt <- read.csv("plotcount.csv")

com_mat<- as.data.frame.matrix(table(dt$Site, dt$Species)) #Sitewise

com_mat2<- as.data.frame.matrix(table(dt$Plotcode, dt$Species)) #Plotwise

write.csv(com_mat, file="sitecommat.csv")
write.csv(com_mat2, file="plotcommat.csv")


# NMDS ------------------------------------------------

#Read community matrix
#abund<-read.csv2("commatrix.csv", sep=",")
abund<-read.csv2("commat2.csv", sep=",")

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

# Anosim test difference between communities ------------------------------

#Adonis analyzes and partitions sums of squares using distance matrices. It can be seen as an ANOVA using distance matrices (analogous to MANOVA - multivariate analysis of variance). Therefore, it is used to test if two or more groups have similar compositions.
BC.dist=vegdist(abund[,3:33], distance="bray")

adonis(BC.dist ~ abund$Landuse, perm=999)

anosim(BC.dist, abund$Landuse, permutations = 1000)

# Nestedness and Turnover -------------------------------------------------

#betapart
dist<-bray.part(abund[,3:33])
bd <- beta.multi.abund(dist[[3]], index.family="bray")
bd
bd <- beta.sample.abund(abund[,3:33], index.family="bray")


abund1 <- abund[1:9,]
abund2 <- abund[10:18,]

dist1 <- beta.pair.abund(abund1[,3:34], index.family = "bray")
dist2 <- beta.pair.abund(abund2[,3:34], index.family = "bray")
mean(dist1[[1]])
mean(dist1[[2]])
mean(dist1[[3]])

mean(dist2[[1]])
mean(dist2[[2]])
mean(dist2[[3]])

core <- beta.sample.abund(abund[,3:34])
dist$beta.ruz.bal

plot(bd)
boxplot(bd)
anova(bd)

dist<-beta.pair.abund(abund[,3:33], index.family="bray")
dist
adonis(dist[[1]]~ abund$Landuse, perm=999)
adonis(dist[[2]]~ abund$Landuse, perm=999)


dist<-bray.part(abund[,3:33])
a<- mean(dist[[1]]) 
b<- mean(dist[[2]])
c<- mean(dist[[3]])




groups <- (abund[,2])

groups <- factor(c(rep(1,9), rep(2,9)), labels = c("commonland", "forest"))

bd<-betadisper(dist[[3]],groups)
bd


plot(bd, main="beta diversity (bray-curtis-index), woody plants")
boxplot(bd)
title(main =" distance of ??-diversity values in relation to their centroids")

anova(bd)

bd

beta.sample.abund(abund[,3:33], index.family="bray", sites=18, samples=100)
oecosimu(abund[,3:33], nestedchecker, "r0")

mod <- betadisper(dist[[3]], groups)
anova(mod)
TukeyHSD(mod)

#adespatial Legendre 2014
install.packages("adespatial")
library(adespatial)
beta.div.comp(abund[,3:34], coef = "S", quant = TRUE, save.abc = FALSE)
?beta.div.comp

adonis(dist[[3]] ~ abund$Landuse, perm=999)
sd(dist[[1]])
mean(dist[[3]])

# Pie chart ---------------------------------------------------------------

data <- data.frame(
  group= c("Turnover", "Nestedness"),
  value=c(0.5026923, 0.09255029)
)

ggplot(data, aes(x="", y=value, fill=group)) +
  geom_bar(stat="identity", width=1, color="white")
+
  coord_polar("y", start=0) +
  
  theme_void() # remove background, grid, numeric labels


#barplot
# data <- data.frame(
#   group= c("Total beta diversity", "Turnover", "Nestedness" ),
#   value=c( 0.5952425, 0.5026923, 0.09255029)
# )

data <- data.frame(
  group= c("Total beta diversity", "Turnover", "Nestedness" ),
  value=c( 0.5836464, 0.518241, 0.06540536)
)

data1 <- data.frame(
  group= c("Total beta diversity", "Turnover", "Nestedness" ),
  value=c( 0.5038292, 0.3437083, 0.1601208)
)

data1 <- data.frame(
  group= c("Turnover", "Nestedness" ),
  value=c( 0.3437083, 0.1601208)
)

data2 <- data.frame(
  group= c("Total beta diversity", "Turnover", "Nestedness" ),
  value=c( 0.5575908, 0.3729472, 0.1846436)
)

data2 <- data.frame(
  group= c("Turnover", "Nestedness" ),
  value=c(0.3729472, 0.1846436)
)

data3 <- data.frame(
  group= c("Turnover", "Nestedness" ),
  value=c(0.7934752, 0.05985987)
)
total <- rbind(data1, data2)

a <- ggplot(data1, aes(x=fct_inorder(group), y=value, fill=group)) +
  geom_bar(stat="identity", width=.5, position = "dodge")+
  theme(axis.text.x=element_blank())+ theme_classic()
a + labs(x = "", y = " ") +
  theme(legend.position="none")

b <- ggplot(data2, aes(x=fct_inorder(group), y=value, fill=group)) +
  geom_bar(stat="identity", position = position_stack())+
  theme(axis.text.x=element_blank())+ theme_classic()
b + labs(x = "", y = "") +
  theme(legend.position="none")

ggplot(total, aes( x= group , y = value, fill = group))+ geom_bar(stat="identity", position = position_stack())


ggplot(ref_by_period_df, aes(x = Hex, y = Qty, fill = Teacher)) + 
  geom_bar(stat="identity") + 
  facet_grid(.~Period) +
  theme(axis.text.x = element_text(angle = 90, size = 8))

library(cowplot)
plot_grid(a, b, labels=c("Open savanna", "Reserve forest"))
?plot_grid


0.7934752


dat <- read.csv("barcomp.csv")

ggplot(dat, aes(Habitat, value)) + 
  geom_bar(aes(fill =div), stat = 'identity', position = "stack", width=0.7) + 
  ylab("Contributions to ?? diversity") +
  theme(axis.text.x=element_blank())+ theme_classic() + 
  theme(legend.title=element_blank()) +scale_fill_brewer(palette="Blues")


a <- ggplot(data3, aes(group, value)) +
  geom_bar(aes(fill= group), stat="identity", position ="stack", width=.5)+ 
  ylab("Contributions to ?? diversity")+
  theme(axis.text.x=element_blank())+ theme_classic()+ theme(legend.title=element_blank()) +scale_fill_brewer(palette="Reds")

a + labs(x = "", y = " ") +
  theme(legend.position="none")

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

write.csv(b, file= "ivig.csv")

data(ifri)
?importancevalue

# Relative abundance ------------------------------------------------------
library(forcats)
relabund <- read.csv("relabund2.csv")
relabund$Species <- as.character(relabund$Species)
relabund$Species<- factor(relabund$Species, levels = relabund$Species)

dodge <- position_dodge(width = 0.8)

p<- ggplot(relabund, aes(x= Species, fill=variable, y= value)) + 
  geom_bar(position = position_dodge(0.8), width = 0.5 , stat="identity")+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))

+scale_colour_manual(name = 'Habitat',
                     values = colors) +
  scale_fill_manual(name = 'Habitat',
                    values = colors) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+ scale_y_continuous(expand = c(0,0)) + geom_text(aes(x=1, y=0.35, label="Stretch it"), vjust=-1) 



p + aes(x = fct_inorder(Species)) + theme(
  panel.background = element_rect(fill = "transparent"), # bg of the panel
  plot.background = element_rect(fill = "transparent", color = NA), # bg of the plot
  # panel.grid.major = element_blank(), # get rid of major grid
  # panel.grid.minor = element_blank(), # get rid of minor grid
  legend.background = element_rect(fill = "transparent")# get rid of legend bg
  # legend.box.background = element_rect(fill = "transparent") # get rid of legend panel bg
  ,axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 9, face="italic"), axis.title.x = element_text(size=10), axis.title.y =element_text(size=10) , legend.title = element_text(size=7), legend.text=element_text(size=7), legend.key.size = unit(0.3, "cm")) + labs(x = "Species", y = "Relative abundance")  + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text= element_text(size = 8)) 


p+ aes(x = fct_inorder(Species)) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), text= element_text(size = 8) ,axis.line = element_line(colour = "black"), axis.text.x = element_text(size = 8, face="italic"), axis.title.x = element_text(size=8), axis.title.y =element_text(size=8)) + labs(x = "Species", y = "Relative abundance") 

ggsave(file="relabundgec.png", width=6, height=4, dpi=300)


# Size structure ----------------------------------------------------------

sizes_g <- read.csv("sizes_cl2.csv")
sizes_f <- read.csv("sizes_rf2.csv")

# ggplot(sizest) + 
#   geom_bar(aes(x = Size.class, fill = factor(Botanical.name)), position = position_dodge(preserve = 'single')) + theme_classic()

par(mfrow=c(2,1))

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

#forest
fit_fw  <- fitdist(Sizes_f, "weibull")
summary(fit_fw)

#grassland
fit_gw  <- fitdist(Sizes_g, "weibull")
summary(fit_gw)

par(mfrow=c(2,1))
plotdist(Sizes_f,demp = TRUE)

denscomp(list(fit_fw), main="Wooded savanna", breaks=15, xlab="Size class in cm")
denscomp(list(fit_gw), main="Open savanna", breaks=15, xlab="Size class in cm")


# Distance decay ----------------------------------------------------------

#Dissimilarity matrix spatial matrix
setwd("C:/Users/vn38fiva/Desktop/chapter 3/chap3analysis/chap3/distdecay")

install.packages("geosphere")
library(geosphere)
df_cl = read.csv("cl.csv", header= TRUE, sep=",")
df_rf = read.csv("rf.csv", header= TRUE, sep=",")

# df = read.csv("2008forestmat.csv", header= TRUE, sep=",")
# df = read.csv("2008nearforestmat.csv", header= TRUE, sep=",")
#abundance data frame
abund = df_cl[,5:ncol(df_cl)]
#abund = df[,4:ncol(df)]

spat.dist<-dist(data.frame(df_cl$Longitude, df_cl$Latitude))

#longitude and latitude 
geo = data.frame(df_cl$Longitude, df_cl$Latitude)


#abundance data frame - bray curtis dissimilarity
dist.abund = vegdist(abund, method = "bray")

#geographic data frame - haversine distance 
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

#abundance vs geographic 
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo          

plot(dist.geo,dist.abund,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
abline(lm(dist.abund ~ dist.geo))

library(magrittr)
lm(dist.abund ~ dist.geo) %>% summary

abund = df_rf[,5:ncol(df_rf)]
#abund = df[,4:ncol(df)]

spat.dist<-dist(data.frame(df_rf$Longitude, df_rf$Latitude))

#longitude and latitude 
geo = data.frame(df_rf$Longitude, df_rf$Latitude)


#abundance data frame - bray curtis dissimilarity
dist.abund = vegdist(abund, method = "bray")

#geographic data frame - haversine distance 
d.geo = distm(geo, fun = distHaversine)
dist.geo = as.dist(d.geo)

#abundance vs geographic 
abund_geo  = mantel(dist.abund, dist.geo, method = "spearman", permutations = 9999, na.rm = TRUE)
abund_geo          

plot(dist.geo,dist.abund,pch=16,cex=0.5,col="black",bty="l",xlab="Spatial distance",ylab="Species composition dissimilarity")
abline(lm(dist.abund ~ dist.geo))

#Subplot nested in site glmer
df_dist <- data.frame(d_abund = dist.abund %>% as.vector,
                      d_m = dist.geo %>% as.vector,
                      site = df$Site)
glmer(dist.abund ~ dist.geo + (1|Site), family=) %>% summary


# SAC ---------------------------------------------------------------------

#Neredicheruvu
library(vegan)
# data(BCI)
# 
# #Vegan data is a SITE x SPECIES matrix
# head(BCI[,1:3])
#Default specaccum
# sp1 <- specaccum(BCI)
# 
# # specaccum w/ method = "collector"
# sp2 <- specaccum(BCI, method = "collector")
# 
# 
# par(mfrow = c(2,2), mar = c(4,2,2,1))
# plot(sp1, main = "Default: plot( specaccum(BCI) )")
# plot(sp2, main = "method = collector")
# plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue", 
#      main = "Default: Prettier CI")

setwd("C:/Users/vn38fiva/Desktop/chapter 3/chap3analysis/chap3/Neredicheruvu")

abund_n<-read.csv2("neredicheruvu.csv", sep=",")
abund_nr <- abund_n[,-1]
sp1 <- specaccum(abund_nr)
sp2 <- specaccum(abund_nr, method = "collector")
par(mfrow = c(2,2), mar = c(4,2,2,1))
# plot(sp1, main = "Default: plot( specaccum(abund_nr) )")
# plot(sp2, main = "method = collector")
plot(sp1, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     main = "Default: Prettier CI")

abund_m<-read.csv2("ng.csv", sep=",")
abund_ng <- abund_m[,5:35]
sp3 <- specaccum(abund_ng)
sp4 <- specaccum(abund_ng, method = "collector")
par(mfrow = c(2,2), mar = c(4,2,2,1))
# plot(sp1, main = "Default: plot( specaccum(abund_nr) )")
# plot(sp2, main = "method = collector")
plot(sp3, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue",
     main = "Default: Prettier CI")
