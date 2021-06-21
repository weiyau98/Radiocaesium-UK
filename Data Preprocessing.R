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


#Import data
setwd('C:/Users/ASUS/Desktop/UM/Research Project/UK data/data')
df <- read.csv('Post_Chernobyl_survey_data_radiocaesium.csv', header=T, na.strings=c('','NA'))
data <- df[df$Sampling_date_yyyymmdd!='1/6/1986',] 
#we will be using data from oct_1986 and mar_1987 only


#Oct 1986
#Cs-137 in vegetation
#Filter plants
plants <- data[data$Classification=='Plants',]

#Drop rows with missing values in Cs-137 in vegetation
plants <- plants %>% drop_na(Cs137_bqkg_DM_veg)

#Drop duplicated rows that have similar values
plants <- plants[plants$code!='BB125G' & plants$code!='N0873G',]

#Arrange row number
rownames(plants) <- 1:nrow(plants)

#Factor variables
plants$Sample <- as.factor(plants$Sample)
plants$Soil_texture_upper <- as.factor(plants$Soil_texture_upper)
plants$Soil_texture_lower <- as.factor(plants$Soil_texture_lower)

#Get variables
date <- c('Sampling_date_yyyymmdd')
loc <- c('Latitude','Longitude')
nvars <- c('Dry_matter_biomass_g','Soil_pH_upper','Depth_upper_soil_layer',
           'LOI_upper_soil_layer','Upper_bulk_density','Soil_pH_lower_soil',
           'Depth_lower_soil')
cvars <- c('Sample_type','Soil_texture_upper','Soil_texture_lower')
conc <- c('Cs137_bqkg_DM_veg')

#Get unique coordinates from latitude and longitude
plants_loc <- select(plants, loc)
plants_loc$coord <- paste(plants_loc$Latitude, plants_loc$Longitude, sep=",")
plants_loc_u <- plants_loc %>% distinct(coord, .keep_all=T)
plants_loc_u <- plants_loc_u[,-3]

#Obtain neighbours list
coordsm <- as.matrix(plants_loc_u)
coordsm <- coordsm[,c('Longitude','Latitude')]
nb <- graph2nb(gabrielneigh(coordsm), sym=TRUE)
g <- s.label(coordsm, nb = nb, pSp.col = "grey", 
             pnb.edge.col = "red", ppoints.cex = 1.5,
             plabels.cex = 0, plot = FALSE)
ADEgS(list(g))

#Obtain edge weight by inverse harvesine distances
invdist <- lapply(nbdists(nb, coordsm, longlat=T), function(x) 1/x)

#Obtain Spatial weight matrix
lw <- nb2listw(nb, glist=invdist, style="W")
can.be.simmed(lw) #Check whether is symmetric

#Moran's eiegenvector map
me <- mem(lw)
s.value(coordsm, me[,4], symbol='circle', ppoint.cex=0.9, porigin.include=F)
dev.copy2pdf(file="mem4.pdf",out.type="cairo", width=10, height=7.5)
s.value(coordsm, me[,137], symbol='circle', ppoint.cex=0.9, porigin.include=F)
dev.copy2pdf(file="mem137.pdf",out.type="cairo", width=10, height=7.5)

#MEM variables selection
#Obtain median for each numerical variable by unique coordinates
plants_nvars <- select(plants, loc, nvars)
plants_nvars <- plants_nvars %>% group_by(Latitude, Longitude) %>% summarise_each(~median(., na.rm=T))
plants_nvars <- join(plants_loc_u, plants_nvars, 
                     by=c('Latitude','Longitude'), type='left') #arrange row according to coordsm

#Impute missing values for data for MEM selection
imputed <- mice(plants_nvars, m=5, maxit=20, seed=123)
plants_nvars_imputed <- complete(imputed, 3)

#Select MEMs
me_sel <- mem.select(plants_nvars_imputed, listw=lw)
me_sel$MEM.select
mem <- as.data.frame(me_sel$MEM.select)

#Final dataset with all variables
coord_mem <- cbind(plants_loc_u, mem)
plants_df <- select(plants, date, loc, nvars, cvars, conc)
plants_df <- join(plants_df, coord_mem, by=c('Latitude','Longitude'), type='left')

#Create indicator function for train(1) and test(0)
plants_df$IND <- ifelse(plants_df$Sampling_date_yyyymmdd=='15/10/1986', '1', '0')
write.csv(plants_df, 'C:/Users/ASUS/Desktop/UM/Research Project/guide/plants.csv', row.names=F)

#Transformation on response variable
#Log transformation
plants_df$Cs137_bqkg_DM_veg <- log(plants_df$Cs137_bqkg_DM_veg, base=10)
write.csv(plants_df, 'C:/Users/ASUS/Desktop/UM/Research Project/guide/plantslog.csv', row.names=F)


#Output text file for GUIDE
setwd('C:/Users/ASUS/Desktop/UM/Research Project/guide')
plants_g <- read.csv('C:/Users/ASUS/Desktop/UM/Research Project/guide/plants.csv')
G <- data.frame(1:ncol(plants_g), colnames(plants_g), 'n')
write.table(G, file='plants.DSC', quote=F, row.names=F)