#Load packages
library(rgdal) 
library(plyr)
library(dplyr)
library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(tidyr)
library(mice)
library(ggplot2)


#Analysis of prediction
#Random forest 
setwd('C:/Users/ASUS/Desktop/UM/Research Project/guide')
plants <- read.table(file='plants.pred', header=T)
colseq = c(rgb(0,114,178,max='255',alpha=150),rgb(230,159,0,max='255',alpha=255))
colnames(plants)
plot(plants[,2:3], ylim=c(0,2500), xlab='Observed', ylab='Predicted', 
     cex=0.75, pch=16, col=colseq[unclass(as.factor(plants$train))], 
     main='Cs-137 in vegetation')
legend(7000, 500, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))
abline(0,1)

plants$res <- plants$observed - plants$predicted
plot(plants$predicted, plants$res, xlab='Predicted Cs-137 in vegetation', ylab='Residuals', 
     cex=0.75, pch=16, col=colseq[unclass(as.factor(plants$train))], 
     main='Residuals vs predicted Cs-137 in vegetation')
legend(250, 6000, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))
abline(0,0)


#Random forest after log transformation on response variable
plantslog <- read.table(file='plantslog.pred', header=T)
colnames(plantslog)
plot(plantslog[,2:3], xlab='Observed', ylab='Predicted', 
     cex=0.75, pch=16, col=colseq[unclass(as.factor(plantslog$train))], 
     main='Cs-137 in vegetation')
legend(0.5, 3, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))
abline(0,1)

plantslog$res <- plantslog$observed - plantslog$predicted
plot(plantslog$predicted, plantslog$res, xlab='Predicted Cs-137 in vegetation', ylab='Residuals', 
     cex=0.75, pch=16, col=colseq[unclass(as.factor(plantslog$train))], 
     main='Residuals vs predicted Cs-137 in vegetation')
legend(1.15, 1.5, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))
abline(0,0)


#Importance scoring of environmental variables
#Original data
par(mar=c(5,10,4,2),las=1)
leg.col <- c("orange","yellow","white")
leg.txt <- c("highly important","likely important","unimportant")
vscore <- read.table("vscore.des",header=TRUE)
score <- vscore$Score
vars <- vscore$Variable
type <- vscore$Type
barcol <- rep("orange",length(vars))
barcol[type == "L"] <- "yellow"
barcol[type == "U"] <- "white"
barplot(rev(score),names.arg=rev(vars),col=rev(barcol),horiz=TRUE,
        xlab="GUIDE importance scores", cex.name=0.55,
        main='Importance scoring of variables')
abline(v=1,col="black",lty=2)
legend("bottomright",legend=leg.txt,fill=leg.col)

#Data after log transformation on response variable
par(mar=c(5,10,4,2),las=1)
leg.col <- c("orange","yellow","white")
leg.txt <- c("highly important","likely important","unimportant")
vscore <- read.table("vscorelog.des",header=TRUE)
score <- vscore$Score
vars <- vscore$Variable
type <- vscore$Type
barcol <- rep("orange",length(vars))
barcol[type == "L"] <- "yellow"
barcol[type == "U"] <- "white"
barplot(rev(score),names.arg=rev(vars),col=rev(barcol),horiz=TRUE,
        xlab="GUIDE importance scores", cex.name=0.55,
        main='Importance scoring of variables')
abline(v=1,col="black",lty=2)
legend("bottomright",legend=leg.txt,fill=leg.col)


#Scatterplot of MEM4 and log CS
setwd('C:/Users/ASUS/Desktop/UM/Research Project/guide')
colseq = c(rgb(0,114,178,max='255',alpha=150),rgb(230,159,0,max='255',alpha=255))
df <- read.csv(file='plantslog.csv')
plot(df$MEM4,df$Cs137_bqkg_DM_veg,main='MEM4 vs log CS',
     xlab='MEM4',ylab='log CS',pch=16,cex=0.75,
     col=colseq[unclass(as.factor(df$IND))])
legend(4, 4, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))

cor.test(df$MEM4,df$Cs137_bqkg_DM_veg,method='pearson')


#Scatterplot of MEM4 and lower soil pH
plot(df$MEM4,jitter(df$Soil_pH_lower_soil, factor=2),main='MEM4 vs lower soil pH',
     xlab='MEM4',ylab='lower soil pH',pch=16,cex=0.75,
     col=colseq[unclass(as.factor(df$IND))])
legend(4, 7, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))

cor.test(df$MEM4,df$Soil_pH_lower_soil,method='pearson')


#Scatterplot of MEM4 and upper soil pH
plot(df$MEM4,jitter(df$Soil_pH_upper, factor=2),main='MEM4 vs upper soil pH',
     xlab='MEM4',ylab='upper soil pH',pch=16,cex=0.75,
     col=colseq[unclass(as.factor(df$IND))])
legend(4, 4, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))

cor.test(df$MEM4,df$Soil_pH_upper,method='pearson')


#Scatterplot of MEM137 and log CS
plot(df$MEM137,df$Cs137_bqkg_DM_veg,main='MEM137 vs log CS',
     xlab='MEM137',ylab='log CS',pch=16,cex=0.75,
     col=colseq[unclass(as.factor(df$IND))])
legend(4, 4, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))

cor.test(df$MEM137,df$Cs137_bqkg_DM_veg,method='pearson')


#Scatterplot of LOI_upper_soil_layer and log CS
df1 <- df[!is.na(df$LOI_upper_soil_layer),]
plot(df1$LOI_upper_soil_layer,df1$Cs137_bqkg_DM_veg,main='Log CS vs LOI_upper_soil_layer',
     xlab='LOI_upper_soil_layer',ylab='log CS (Bq/kg)',pch=16,cex=0.75,
     col=colseq[unclass(as.factor(df1$IND))])
legend(80, 1, pch=16, pt.cex=1.2, pt.lwd=2,
       col=colseq, c('test','train'))

cor.test(df1$LOI_upper_soil_layer,df1$Cs137_bqkg_DM_veg,method='pearson')


#Residuals visualization
#import shapefile
setwd('C:/Users/ASUS/Desktop/UM/Research Project/UK data/data')
shapefile <- readOGR(dsn='C:/Users/ASUS/Desktop/UM/Research Project/gadm36_GBR_shp',
                     layer='gadm36_GBR_2')

#fortify shapefile
mapdata <- fortify(shapefile)

dfres <- cbind(df,plantslog$res)
colnames(dfres)[which(names(dfres) == 'plantslog$res')] <- 'res'
mybreaks <- c(-1.6,-1.2,-0.8,-0.4,0.4,0.8,1.2,1.6)
rdbu <- c('#b2182b','#d6604d','#f4a582','#fddbc7','#d1e5f0','#92c5de','#4393c3','#2166ac')

ggplot(dfres, aes(x=Longitude, y=Latitude))+ 
  geom_point(aes(colour=res), cex=2)+
  scale_colour_gradientn(name='Normalised residuals', colours=rdbu, breaks=mybreaks, guide='colourbar', aesthetics='colour')+
  geom_path(data=mapdata,aes(x=long, y=lat,group=group), colour="grey70")+
  xlim(-10,+2.5)+
  coord_fixed()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())


#Boxplot for Cs137_bqkg_DM_veg and Sample_type
ggplot(aes(x=Sample_type, y=Cs137_bqkg_DM_veg, fill=Sample_type), data=df)+
  geom_boxplot(alpha=0.8)+
  scale_fill_manual(values = c('green3','orange','dodgerblue3'))+
  geom_jitter(color="black", pch=1, size=0.9, alpha=0.8, width=0.1)+
  labs(y='log Cs', x='Sample type')


#ANOVA on logCs and Sample_type
anova <- aov(Cs137_bqkg_DM_veg~Sample_type, data=df)
summary(anova)		#p-value<0.05


#ANOVA on MEM4 and Sample_type
anova <- aov(MEM4~Sample_type, data=df)
summary(anova)		#p-value>0.05


#Distribution of plants 
ggplot(aes(x=Longitude, y=Latitude, color=Sample_type), data=df)+
  geom_point(alpha=0.7, cex=2)+
  scale_colour_manual(values = c('green3','orange','dodgerblue3'), aesthetics='colour')+
  labs(colour='Sample type')


#Scatterplot of upper soil pH and concentration of Cs-137 in vegetation
plot(jitter(df$Soil_pH_upper), df$Cs137_bqkg_DM_veg, xlab="Upper soil pH", 
     ylab="log Cs", xlim=c(2,10), cex=0.8, pch=c(15,16,17)[unclass(as.factor(df$Sample_type))], 
     col=c('green3','orange','dodgerblue3')[unclass(as.factor(df$Sample_type))])
legend(8, 3.75, pch=c(15,16,17), pt.cex=1.2, pt.lwd=2,
       col=c('green3','orange','dodgerblue3'), c('Bracken','Grass','Heather'))

cor.test(df$Soil_pH_upper,df$Cs137_bqkg_DM_veg,method='pearson')


#Scatterplot of lower soil pH and concentration of Cs-137 in vegetation
plot(jitter(df$Soil_pH_lower_soil, factor=6), df$Cs137_bqkg_DM_veg, xlab="Lower soil pH", 
     ylab="log Cs", xlim=c(2,10), cex=0.8, pch=c(15,16,17)[unclass(as.factor(df$Sample_type))],
     col=c('green3','orange','dodgerblue3')[unclass(as.factor(df$Sample_type))])
legend(8, 3.75, pch=c(15,16,17), pt.cex=1.2, pt.lwd=2,
       col=c('green3','orange','dodgerblue3'), c('Bracken','Grass','Heather'))

cor.test(df$Soil_pH_lower_soil,df$Cs137_bqkg_DM_veg,method='pearson')






