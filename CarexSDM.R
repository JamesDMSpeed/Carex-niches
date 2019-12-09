#SDM for Carex lep and Carex jem

library(raster)
library(rasterVis)
library(rgdal)
library(sp)
library(sdm)
library(gridExtra)


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

#Convert temp to degrees
predvars$bio10_16<-predvars$bio10_16/10
#Convert soil pH to decimal ph
predvars$SoilpH<-predvars$SoilpH/10


norwayP<-spTransform(norway,CRS=crs(predvars))


#Importing Species data
carlep<-read.csv('Species data/C-lep-gbif-final.csv',header=T,sep=';')
View(carlep)
carjem<-read.csv('Species data/C-jemt-gbif-final.csv',header=T,sep=';')
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

#Single object
carex <- rbind(carlep_utm[,names(carlep_utm)%in%c("decimalLatitude","decimalLongitude")],carjem_utm[,names(carjem_utm)%in%c("decimalLatitude","decimalLongitude")])
carex$species<-c(rep('Carex_lepidocarpa',times=nrow(carlep_utm)),rep('Carex_jemtlandica',times=nrow(carjem_utm)))


pt<-levelplot(predvars$bio10_16,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main=expression('Mean summer temperature'~(degree~C)),
          colorkey=list(space='bottom',title=''))+
  layer(sp.polygons(norwayP,col=grey(0.5)))
 #layer(sp.points(carlep_utm,col=1,pch=16,cex=0.5))

pps<-levelplot(predvars$bio15_16,scales=list(draw=F),margin=F,par.settings='RdBuTheme',main='Precipitation seasonality',
              colorkey=list(space='bottom',title=expression('')))+
  layer(sp.polygons(norwayP,col=grey(0.5)))

pph<-levelplot(mask(predvars$SoilpH,predvars$bio10_16),scales=list(draw=F),margin=F,par.settings='RdBuTheme',main='Soil pH',
               colorkey=list(space='bottom',title=expression('')))+
  layer(sp.polygons(norwayP,col=grey(0.5)))


pp<-levelplot(predvars$bio12_16,scales=list(draw=F),margin=F,par.settings='RdBuTheme',main=expression('Mean annual precipitation'~(mm)),
                colorkey=list(space='bottom'))+
  layer(sp.polygons(norwayP,col=grey(0.5)))
  
pL<-levelplot(predvars[[1]]*0,col.regions=grey(1),colorkey=F,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main='Carex lepidocarpa')+
  layer(sp.polygons(norwayP,col=grey(0.5)))+
  layer(sp.points(carlep_utm,col=1,pch=16,cex=0.5))

pJ<-levelplot(predvars[[1]]*0,col.regions=grey(1),colorkey=F,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main='Carex jemtlandica')+
  layer(sp.polygons(norwayP,col=grey(0.5)))+
  layer(sp.points(carjem_utm,col=1,pch=16,cex=0.5))


#Background data
background_dat<-read.table('Species data/BackgroundBiasCorrected.txt',header=T)

background_utm<-SpatialPoints(background_dat,proj4string=CRS('+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0'))

pB<-levelplot(predvars[[1]]*0,col.regions=grey(1),colorkey=F,scales=list(draw=F),margin=F,par.settings='YlOrRdTheme',main='Background data')+
  layer(sp.polygons(norwayP,col=grey(0.5)))+
  layer(sp.points(background_utm,pch=1,cex=0.1,col=1))

tiff('Figures/AllData.tif',width = 1600,height = 1000,pointsize = 20,units = 'px')
grid.arrange(pt,pp,pps,pph,pL,pJ,pB,ncol=4)
dev.off()


#SDM
#Make the sdm dataset
sdmdataset<-sdmData(species~
                    +bio10_16+bio12_16+bio15_16+SoilpH,
                    train=carex,predictors=predvars,bg=list(sample(background_utm,1000),remove=T))
sdmdataset

#SDM
sdm_carex<-sdm(Carex_lepidocarpa+Carex_jemtlandica~
                            bio10_16+bio12_16+bio15_16+SoilpH,
                            data=sdmdataset,
                            methods=c('glm','gam','rf','gbm','mda','fda','brt'),
                            replication=c('cv'),cv.folds=5)
sdm_carex

modevalcarex<-cbind(sdm_carex@run.info,getEvaluation(sdm_carex))
with(modevalcarex,tapply(AUC,species,mean))
with(modevalcarex,tapply(AUC,species,sd))

#Variables importances
varimplist<-list()
for (i in 1:max(sdm_carex@run.info$modelID)){
  ifelse(sdm_carex@run.info$success[i]==TRUE,
         {
         ifelse(!is.null(getVarImp(sdm_carex,id=i)),
                  varimplist[[i]]<-getVarImp(sdm_carex,id=i)@varImportance,
                  varimplist[[i]]<-df1)
           varimplist[[i]]$species<-sdm_carex@run.info$species[i]
           varimplist[[i]]$method<-sdm_carex@run.info$method[i]
           varimplist[[i]]$repid<-sdm_carex@run.info$replicationID[i]}
         ,print(paste('Model failiure run ',i)))
}
CarexVarImp<-do.call('rbind',varimplist)

#Summarise by species
sem<-function(x)sd(x,na.rm=T)/sqrt(length(!is.na(x)))
varimpsp_mean<-with(CarexVarImp,tapply(corTest,list(species,variables),mean))  
varimpsp_sem<-with(CarexVarImp,tapply(corTest,list(species,variables),sem))  

tiff('Figures/VarImp.tif',width=800,height=600,pointsize = 20)
par(mar=c(5,5,1,1))
b1<-barplot(varimpsp_mean,beside=T,ylab='Variable importance',names.arg=c('MST','MAP','Precipitation \n seasonality','Soil pH'),las=1,ylim=c(0,0.80),
            legend.text=sub("_"," ",rownames(varimpsp_mean)))     
arrows(b1,varimpsp_mean+varimpsp_sem,b1,varimpsp_mean-varimpsp_sem,length=0.05,code=3,angle=90)
dev.off()

#Response curves
responsecurvelist<-list()
for (i in 1:length(levels(as.factor(sdm_carex@run.info$species)))){
  responsecurvelist[[i]]<-getResponseCurve(sdm_carex,id=sdm_carex@run.info$modelID[sdm_carex@run.info$species==levels(as.factor(sdm_carex@run.info$species))[i]
                                                                                     &sdm_carex@run.info$method%in% c('glm','gam','brt')]
                                           ,mean=T,main=levels(as.factor(sdm_carex@run.info$species))[i])
}

tiff('Figures/ResponseCurves.tif',width = 1200,height=1200,units='px',pointsize = 20)
par(mfrow=c(2,2))
par(mar=c(5,5,1,1))
plot(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
       ,type='l',main=levels(carex$species)[i],xlab=expression('MST'~(degree~C)),ylab='Response',las=1,ylim=c(0,0.6))
lines(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
        +apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio10_16[,1],apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,mean)
        -apply(responsecurvelist[[1]]@response$bio10_16[,2:ncol(responsecurvelist[[1]]@response$bio10_16)],1,sd),lty=2)
legend('topl',lty=1,col=c(1,2),c('Carex lepidocarpa','Carex jemtlandica'))
lines(responsecurvelist[[2]]@response$bio10_16[,1],apply(responsecurvelist[[2]]@response$bio10_16[,2:ncol(responsecurvelist[[2]]@response$bio10_16)],1,mean)
     ,type='l',main=levels(carex$species)[i],col=2)
lines(responsecurvelist[[2]]@response$bio10_16[,1],apply(responsecurvelist[[2]]@response$bio10_16[,2:ncol(responsecurvelist[[2]]@response$bio10_16)],1,mean)
      +apply(responsecurvelist[[2]]@response$bio10_16[,2:ncol(responsecurvelist[[2]]@response$bio10_16)],1,sd),lty=2,col=2)
lines(responsecurvelist[[2]]@response$bio10_16[,1],apply(responsecurvelist[[2]]@response$bio10_16[,2:ncol(responsecurvelist[[2]]@response$bio10_16)],1,mean)
      -apply(responsecurvelist[[2]]@response$bio10_16[,2:ncol(responsecurvelist[[2]]@response$bio10_16)],1,sd),lty=2,col=2)

plot(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
     ,type='l',main=levels(carex$species)[i],xlab='MAP (mm)',ylab='Response',las=1,ylim=c(0,0.6))
lines(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
      +apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio12_16[,1],apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,mean)
      -apply(responsecurvelist[[1]]@response$bio12_16[,2:ncol(responsecurvelist[[1]]@response$bio12_16)],1,sd),lty=2)

lines(responsecurvelist[[2]]@response$bio12_16[,1],apply(responsecurvelist[[2]]@response$bio12_16[,2:ncol(responsecurvelist[[2]]@response$bio12_16)],1,mean)
      ,type='l',main=levels(carex$species)[i],col=2)
lines(responsecurvelist[[2]]@response$bio12_16[,1],apply(responsecurvelist[[2]]@response$bio12_16[,2:ncol(responsecurvelist[[2]]@response$bio12_16)],1,mean)
      +apply(responsecurvelist[[2]]@response$bio12_16[,2:ncol(responsecurvelist[[2]]@response$bio12_16)],1,sd),lty=2,col=2)
lines(responsecurvelist[[2]]@response$bio12_16[,1],apply(responsecurvelist[[2]]@response$bio12_16[,2:ncol(responsecurvelist[[2]]@response$bio12_16)],1,mean)
      -apply(responsecurvelist[[2]]@response$bio12_16[,2:ncol(responsecurvelist[[2]]@response$bio12_16)],1,sd),lty=2,col=2)

plot(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
     ,type='l',main=levels(carex$species)[i],xlab='Precipitation seasonality',ylab='Response',las=1,ylim=c(0,0.6))
lines(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
      +apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$bio15_16[,1],apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,mean)
      -apply(responsecurvelist[[1]]@response$bio15_16[,2:ncol(responsecurvelist[[1]]@response$bio15_16)],1,sd),lty=2)

lines(responsecurvelist[[2]]@response$bio15_16[,1],apply(responsecurvelist[[2]]@response$bio15_16[,2:ncol(responsecurvelist[[2]]@response$bio15_16)],1,mean)
      ,type='l',main=levels(carex$species)[i],col=2)
lines(responsecurvelist[[2]]@response$bio15_16[,1],apply(responsecurvelist[[2]]@response$bio15_16[,2:ncol(responsecurvelist[[2]]@response$bio15_16)],1,mean)
      +apply(responsecurvelist[[2]]@response$bio15_16[,2:ncol(responsecurvelist[[2]]@response$bio15_16)],1,sd),lty=2,col=2)
lines(responsecurvelist[[2]]@response$bio15_16[,1],apply(responsecurvelist[[2]]@response$bio15_16[,2:ncol(responsecurvelist[[2]]@response$bio15_16)],1,mean)
      -apply(responsecurvelist[[2]]@response$bio15_16[,2:ncol(responsecurvelist[[2]]@response$bio15_16)],1,sd),lty=2,col=2)

plot(responsecurvelist[[1]]@response$SoilpH[,1],apply(responsecurvelist[[1]]@response$SoilpH[,2:ncol(responsecurvelist[[1]]@response$SoilpH)],1,mean)
     ,type='l',main=levels(carex$species)[i],xlab='Soil pH',ylab='Response',las=1,ylim=c(0,0.6))
lines(responsecurvelist[[1]]@response$SoilpH[,1],apply(responsecurvelist[[1]]@response$SoilpH[,2:ncol(responsecurvelist[[1]]@response$SoilpH)],1,mean)
      +apply(responsecurvelist[[1]]@response$SoilpH[,2:ncol(responsecurvelist[[1]]@response$SoilpH)],1,sd),lty=2)
lines(responsecurvelist[[1]]@response$SoilpH[,1],apply(responsecurvelist[[1]]@response$SoilpH[,2:ncol(responsecurvelist[[1]]@response$SoilpH)],1,mean)
      -apply(responsecurvelist[[1]]@response$SoilpH[,2:ncol(responsecurvelist[[1]]@response$SoilpH)],1,sd),lty=2)

lines(responsecurvelist[[2]]@response$SoilpH[,1],apply(responsecurvelist[[2]]@response$SoilpH[,2:ncol(responsecurvelist[[2]]@response$SoilpH)],1,mean)
      ,type='l',main=levels(carex$species)[i],col=2)
lines(responsecurvelist[[2]]@response$SoilpH[,1],apply(responsecurvelist[[2]]@response$SoilpH[,2:ncol(responsecurvelist[[2]]@response$SoilpH)],1,mean)
      +apply(responsecurvelist[[2]]@response$SoilpH[,2:ncol(responsecurvelist[[2]]@response$SoilpH)],1,sd),lty=2,col=2)
lines(responsecurvelist[[2]]@response$SoilpH[,1],apply(responsecurvelist[[2]]@response$SoilpH[,2:ncol(responsecurvelist[[2]]@response$SoilpH)],1,mean)
      -apply(responsecurvelist[[2]]@response$SoilpH[,2:ncol(responsecurvelist[[2]]@response$SoilpH)],1,sd),lty=2,col=2)
dev.off()

#Predictions
leppred<-predict(sdm_carex,predvars,filename='ModelPredictions/carlep',species='Carex_lepidocarpa',mean=T)
leppred  
jempred<-predict(sdm_carex,predvars,filename='ModelPredictions/carjem',species='Carex_jemtlandica',mean=T)
jempred  
leppred_mean<-calc(leppred,mean)
jempred_mean<-calc(jempred,mean)
levelplot(stack(leppred_mean,jempred_mean),scales=list(draw=F),names.attr=c('C. lepidocarpa','C. jemtlandica'))

#Ensembling
ens_lep<-ensemble(sdm_carex,newdata = predvars,filename = 'ModelPredictions/lep_ensemble',
                  setting=list(method='weighted',stat='AUC',
                               id=sdm_carex@run.info$modelID[sdm_carex@run.info$species=='Carex_lepidocarpa']),overwrite=TRUE)
ens_jem<-ensemble(sdm_carex,newdata = predvars,filename = 'ModelPredictions/jem_ensemble',
                  setting=list(method='weighted',stat='AUC',
                               id=sdm_carex@run.info$modelID[sdm_carex@run.info$species=='Carex_jemtlandica']),overwrite=T)

tiff('Figures/HSMMap.tif',width=1200,height=800,units='px',pointsize = 20,res=300)
levelplot(stack(ens_lep,ens_jem),scales=list(draw=F),names.attr=c('C. lepidocarpa','C. jemtlandica'),par.settings='YlOrRdTheme',
          cex=1.5)+
  layer(sp.polygons(norwayP,col=grey(0.5)))
dev.off()

#Niche (MST and precip season)
niche(predvars,ens_lep,c('bio10_16','bio15_16'))
nl<-niche(predvars,ens_lep,c('bio10_16','bio15_16'),plot=F)
nj<-niche(predvars,ens_jem,c('bio10_16','bio15_16'),plot=F)

#Stack up niches
s1<-stack(nl@nicheRaster$niche,nj@nicheRaster$niche)
levelplot(s1,par.settings='YlOrRdTheme')
#Set extent as actual climate variables
extent(s1)<-stack(nl@scaleParams)[,1]
tiff('Figures/Niches.tif',width=800,height = 800,res=150)
levelplot(s1,par.settings='YlOrRdTheme',
          xlab=expression('MST'~(degree~C)),ylab='Precipitation seasonality',
          names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)
dev.off()

niche(predvars,ens_lep,c('bio10_16','bio12_16'))
npl<-niche(predvars,ens_lep,c('bio10_16','bio12_16'),plot=F)
npj<-niche(predvars,ens_jem,c('bio10_16','bio12_16'),plot=F)
sp<-stack(npl@nicheRaster$niche,npj@nicheRaster$niche)
levelplot(sp,par.settings='YlOrRdTheme')
#Set extent as actual climate variables
extent(sp)<-stack(npl@scaleParams)[,1]
levelplot(sp,par.settings='YlOrRdTheme',
          xlab=expression('MST'~(degree~C)),ylab='Precipitation (mm)',
          names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)

lp1<-levelplot(sp,par.settings='YlOrRdTheme',
               xlab=expression('MST'~(degree~C)),ylab='Precipitation (mm)',
               names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)
lp2<-levelplot(s1,par.settings='YlOrRdTheme',
               xlab=expression('MST'~(degree~C)),ylab='Precipitation seasonality',
               names.attr=c('C. lepidocarpa','C. jemtlandica'),cex=0.8)
#tiff('Figures/Niches2.tif',height = 1200,width = 1000,units='px',res=100)
grid.arrange(lp1,lp2)
#dev.off()
