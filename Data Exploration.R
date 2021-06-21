#Load packages
library(rgdal) 
library(plyr)
library(dplyr)
library(ade4)
library(adespatial)
library(adegraphics)
library(spdep)
library(ggplot2)


#Import data
setwd('C:/Users/ASUS/Desktop/UM/Research Project/UK data/data')
df <- read.csv('Post_Chernobyl_survey_data_radiocaesium.csv', header=T, na.strings=c('','NA'))
june_1986 <- df[df$Sampling_date_yyyymmdd=='1/6/1986',] 
oct_1986 <- df[df$Sampling_date_yyyymmdd=='15/10/1986',] 
mar_1987 <- df[df$Sampling_date_yyyymmdd=='1/3/1987',] 


#Import shapefile
shapefile <- readOGR(dsn='C:/Users/ASUS/Desktop/UM/Research Project/gadm36_GBR_shp',
                     layer='gadm36_GBR_2')


#Fortify shapefile
mapdata <- fortify(shapefile)


#Plot map with points
#Cs137_bqkg_DM_veg 
#june_1986 #Grass
ggplot(june_1986, aes(x=Longitude, y=Latitude))+ 
  ggtitle("Concentration of Cs-137 in dry matter vegetation on June 1986")+
  geom_point(aes(colour=Cs137_bqkg_DM_veg), cex=2)+
  scale_colour_gradient(limits=c(0,20000),low = "orange", high = "black", guide = "colourbar", aesthetics = "colour", na.value="deepskyblue", name="Cs-137(Bq/kg)")+
  geom_path(data=mapdata,aes(x=long, y=lat,group=group), colour="grey70")+
  xlim(-10,+2.5)+
  coord_fixed()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(axis.title.x=element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank())+
  theme(axis.title.y=element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank())

#Cs-137 concentration in vegetation from june 1986 to mar 1987
plants <- df
plants <- plants[plants$Sample_type=='Grass',c('Sampling_date_yyyymmdd',"Latitude","Longitude","Cs137_bqkg_DM_veg")]
plants <- plants %>%
  mutate(Sampling_date_yyyymmdd=case_when(Sampling_date_yyyymmdd == '1/3/1987' ~ 'March 1987',
                                          Sampling_date_yyyymmdd == '1/6/1986' ~ 'June 1986',
                                          Sampling_date_yyyymmdd == '15/10/1986' ~ 'October 1986'))
plants$Sampling_date_yyyymmdd <- factor(plants$Sampling_date_yyyymmdd, levels = c("June 1986", "October 1986", "March 1987"))
levels(plants$Sampling_date_yyyymmdd)

ggplot(aes(x=Longitude, y=Latitude, size=Cs137_bqkg_DM_veg, color=Sampling_date_yyyymmdd), data=plants)+
  geom_point(alpha=0.65)+
  scale_size(range=c(2,10), name="Cs-137 (Bq/kg)")+
  scale_colour_manual(values = c("#1f78b4","#33a02c","#b2df8a"), aesthetics='colour')+
  labs(colour='Sampling date')
