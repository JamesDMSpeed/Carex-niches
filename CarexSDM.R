#SDM for Carex lep and Carex jem

library(raster)
library(rasterVis)
library(rgdal)
library(sp)
library(sdm)


#Norway outline
norway<-getData('GADM',level=0,country='NOR')

#Environmental variables
list.files('Data/')

#Stack at 1km resolution
predvars<-stack(list.files('Data/',full.names=T))
predvars
#Tidy layer names
names(predvars)<-sub("X_","",names(predvars))
plot(predvars)

norwayP<-spTransform(norway,CRS=crs(predvars))


#Importing Species data
carlep<-read.csv('Species data/C-lep-gbif-filtered.csv',header=T,sep=';')
View(carlep)
carjem<-read.csv('Species data/C-jemt-gbif-filtered.csv',header=T,sep=';')
View(carjem)
str(carlep)
str(carjem)
dim(carlep)
dim(carjem)

#Convert to spatial points dataframe
carlep_sp<-SpatialPointsDataFrame(cbind(carlep$decimalLongitude,carlep$decimalLatitude),carlep,proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
carjem_sp<-SpatialPointsDataFrame(cbind(carjem$decimalLongitude,carjem$decimalLatitude),carjem,proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))
#Project to utm grid
carlep_utm<-spTransform(carlep_sp,CRS=crs(predvars))
carjem_utm<-spTransform(carjem_sp,CRS=crs(predvars))


levelplot(predvars[[1]],scales=list(draw=F),margin=F,par.settings='YlOrRdTheme')+
  layer(sp.polygons(norwayP,col=grey(0.5)))+
  layer(sp.points(carlep_utm))+
  layer(sp.points(carjem_utm,col=1,pch=1))
  