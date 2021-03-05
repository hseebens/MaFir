
graphics.off()
rm(list=ls())

library(sf)
library(shape)

setwd("/home/hanno/Bioinvasion/InvAccu/MarineExpansion/MaFir_Workflow/")

source("../../../SpatialSpread/owncolleg.r") # adjusted color legend


### load data ###############################################################

# all_regspec_fr <- read.table("Data/FirstRecords_TerrMarRegions.csv",sep=";",header=T,stringsAsFactors = F)
all_regspec_fr <- read.table("Data/FirstRecords_TerrMarRegions_min3.csv",sep=";",header=T,stringsAsFactors = F)
# all_regspec_fr <- read.table("../../Others/EkinInternship/FirstRecords_TerrMarRegions_min5.csv",sep=";",header=T,stringsAsFactors = F)

## Polygon file of marine and terrestrial regions
regions <- st_read(dsn="Data/Input/Shapefiles",layer="terraqua",stringsAsFactors = F)
marine_regions <- st_read(dsn="Data/Input/Shapefiles",layer="MEOW_NoOverlap_combined",stringsAsFactors = F)

## change projection to Robinson
crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m" # Robinson projection
regions <- st_wrap_dateline(regions, options = c("WRAPDATELINE=YES"))
regions <- st_transform(regions,crs)
marine_regions <- st_wrap_dateline(marine_regions, options = c("WRAPDATELINE=YES"))
marine_regions <- st_transform(marine_regions,crs)

regions$Realm <- NA
regions$Realm[!is.na(regions$ECOREGION)] <- "marine"
regions$Realm[!is.na(regions$Region)] <- "terrestrial"

regions$Region[!is.na(regions$ECOREGION)] <- regions$ECOREGION[!is.na(regions$ECOREGION)]
# regions$Region[!is.na(regions$featurecla)] <- regions$name[!is.na(regions$featurecla)]
regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)




### plot results #######################################################
terr_timeseries <- as.data.frame(table(subset(all_regspec_fr,Realm=="terrestrial")$FirstRecord),stringsAsFactors = F)
marine_timeseries <- as.data.frame(table(subset(all_regspec_fr,Realm=="marine")$FirstRecord),stringsAsFactors = F)
fresh_timeseries <- as.data.frame(table(subset(all_regspec_fr,Realm=="freshwater")$FirstRecord),stringsAsFactors = F)

x11(width=10,height=3)
op <- par(oma=c(1,1,0,0),mfrow=c(1,3),mar=c(2,2,2,1),mgp=c(1.5,0.5,0),cex=1,las=1,tck=-0.03)

plot(terr_timeseries,xlim=c(1500,2010),main="Terrestrial",ylab="",xlab="")
abline(v=2010,col="gray")
lines(terr_timeseries[,1],runmed(terr_timeseries[,2],k=13),col="red",lwd=2)
mtext("Alien species numbers",side=2,las=0,line=2)
mtext("Years",side=1,las=0,line=1.5)

plot(marine_timeseries,xlim=c(1500,2010),main="Marine",ylab="",xlab="")
abline(v=2010,col="gray")
lines(marine_timeseries[,1],runmed(marine_timeseries[,2],k=13),col="red",lwd=2)
mtext("Years",side=1,las=0,line=1.5)

plot(fresh_timeseries,xlim=c(1500,2010),main="Freshwater",ylab="",xlab="")
abline(v=2010,col="gray")
lines(fresh_timeseries[,1],runmed(fresh_timeseries[,2],k=13),col="red",lwd=2)
mtext("Years",side=1,las=0,line=1.5)
par(op)




x11(width=4,height=3)
op <- par(mar=c(3,3,2,1),mgp=c(1.5,0.5,0),cex=1,las=1,tck=-0.03)

plot(terr_timeseries[,1],cumsum(terr_timeseries[,2])/max(cumsum(terr_timeseries[,2])),xlim=c(1500,2010),ylab="",xlab="",lty=1,type="l",lwd=2)
lines(marine_timeseries[,1],cumsum(marine_timeseries[,2])/max(cumsum(marine_timeseries[,2])),lwd=2,col="blue")
lines(fresh_timeseries[,1],cumsum(fresh_timeseries[,2])/max(cumsum(fresh_timeseries[,2])),lwd=2,col="red")

mtext("Years",side=1,las=0,line=1.5)
par(op)


## world map ###################################
library(RColorBrewer)

nspec_reg <- aggregate(speciesKey ~ Region + Realm,data=all_regspec_fr,FUN=length)
colnames(nspec_reg)[dim(nspec_reg)[2]] <- "nSpec"
# nspec_reg$nSpec[nspec_reg$nSpec>2000] <- 2000


spatial_nspec <- merge(regions,nspec_reg,by=c("Region","Realm"),all.x=T)
# spatial_nspec$nSpec[is.na(spatial_nspec$nSpec)] <- 0
# spatial_nspec <- spatial_nspec[!is.na(spatial_nspec$x),]


### colour coding ######
spatial_nspec$nSpec <- log(spatial_nspec$nSpec+1)

ind_terr <- spatial_nspec$Realm=="terrestrial" | spatial_nspec$Realm=="freshwater" 
spatial_nspec$col_norm[ind_terr] <- round(spatial_nspec$nSpec[ind_terr]/max(spatial_nspec$nSpec[ind_terr],na.rm=T)*99)+1
cols <- colorRampPalette((brewer.pal(n=9,name="YlOrRd")[1:9]))(max(spatial_nspec$col_norm[ind_terr],na.rm=T))
# cols <- rev(colorRampPalette(c("brown4","brown3","orange","yellow","darkolivegreen3","lightblue"))(max(data_regs$col_norm,na.rm=T)))
# cols <- colorRampPalette(c("yellow","yellow","yellow","orange","darkred"))(max(n_sp$Freq_norm)) # 
# cols <- (colorRampPalette(cols)(max(n_sp$Freq_norm)))

spatial_nspec$col <- NA
spatial_nspec[ind_terr,]$col <- cols[spatial_nspec[ind_terr,]$col_norm]
spatial_nspec[ind_terr,]$col[is.na(spatial_nspec[ind_terr,]$col_norm)] <- grey(0.95)

ind_marine <- spatial_nspec$Realm=="marine"
# spatial_nspec$nSpec[ind_marine][spatial_nspec$nSpec[ind_marine]>150] <- 150
spatial_nspec$col_norm[ind_marine] <- round(spatial_nspec$nSpec[ind_marine]/max(spatial_nspec$nSpec[ind_marine],na.rm=T)*99)+1
cols <- colorRampPalette((brewer.pal(n=9,name="YlGnBu")[1:9]))(max(spatial_nspec$col_norm[ind_marine],na.rm=T))

spatial_nspec[ind_marine,]$col <- cols[spatial_nspec[ind_marine,]$col_norm]
spatial_nspec[ind_marine,]$col[is.na(spatial_nspec[ind_marine,]$col_norm)] <- grey(0.95)

# spatial_nspec <- spatial_nspec[order(spatial_nspec$x),]
# cbind(regions$Region,spatial_nspec$Region)


x11(width=8,height=3.3)
# png("../Figures/Worldmap_NumberTaxa.png",unit="in",width=8,height=3.3,res=300)
# png("../Figures/Worldmap_NumberTaxa_GRIISInvasive.png",unit="in",width=8,height=3.3,res=300)
layout(matrix(1:2,nc=2),widths=c(0.85,0.15))
op <- par(mar=c(0,0,0,0),las=1,cex=0.9,tck=-0.02,mgp=c(2,0.3,0))

plot(st_geometry(spatial_nspec),col=spatial_nspec$col,lwd=0.5)

plot(1:10,axes=F,type="n")
owncolleg(posx=c(0.1,0.2),posy=c(0.1,0.95),col=cols,cex=1.2,digit=0,zlim=c(1,4000))#,">2000"
text(x=0.8,y=0.6,labels="Number of alien taxa",srt=270,cex=1.2)
# text(x=0.8,y=0.6,labels="Number of invasive species",srt=270,cex=1.2)
mtext(">",side=3,line=-1.27,adj=0.27)
par(op)
# dev.off()



### Plot marine regions ################################################################

# nspec_reg <- aggregate(speciesKey ~ Region + Realm,data=all_regspec_fr,FUN=length)
# colnames(nspec_reg)[dim(nspec_reg)[2]] <- "nSpec"
# # nspec_reg$nSpec[nspec_reg$nSpec>2000] <- 2000

nspec_reg_mar <- subset(nspec_reg,Realm=="marine")

spatial_nspec <- merge(marine_regions,nspec_reg_mar,by.x="ECOREGION",by.y="Region",all.x=T)
# spatial_nspec$nSpec[is.na(spatial_nspec$nSpec)] <- 0
# spatial_nspec <- spatial_nspec[!is.na(spatial_nspec$x),]

### colour coding ######
spatial_nspec$nSpec <- log10(spatial_nspec$nSpec+1)

# ind_terr <- spatial_nspec$Realm=="terrestrial" | spatial_nspec$Realm=="freshwater" 
# spatial_nspec$col_norm[ind_terr] <- round(spatial_nspec$nSpec[ind_terr]/max(spatial_nspec$nSpec[ind_terr],na.rm=T)*99)+1
# cols <- colorRampPalette((brewer.pal(n=9,name="YlOrRd")[1:9]))(max(spatial_nspec$col_norm[ind_terr],na.rm=T))
# # cols <- rev(colorRampPalette(c("brown4","brown3","orange","yellow","darkolivegreen3","lightblue"))(max(data_regs$col_norm,na.rm=T)))
# # cols <- colorRampPalette(c("yellow","yellow","yellow","orange","darkred"))(max(n_sp$Freq_norm)) # 
# # cols <- (colorRampPalette(cols)(max(n_sp$Freq_norm)))
# 
# spatial_nspec$col <- NA
# spatial_nspec[ind_terr,]$col <- "white" # cols[spatial_nspec[ind_terr,]$col_norm]
# spatial_nspec[ind_terr,]$col[is.na(spatial_nspec[ind_terr,]$col_norm)] <- "grey"#grey(0.95)

# ind_marine <- spatial_nspec$Realm=="marine"
# spatial_nspec$nSpec[ind_marine][spatial_nspec$nSpec[ind_marine]>150] <- 150
spatial_nspec$col_norm <- round(spatial_nspec$nSpec/max(spatial_nspec$nSpec,na.rm=T)*99)+1
cols <- colorRampPalette((brewer.pal(n=9,name="YlGnBu")[3:9]))(max(spatial_nspec$col_norm,na.rm=T))

spatial_nspec$col <- cols[spatial_nspec$col_norm]
spatial_nspec$col[is.na(spatial_nspec$col_norm)] <- "white" # grey(0.95)

# spatial_nspec <- spatial_nspec[order(spatial_nspec$x),]
# cbind(regions$Region,spatial_nspec$Region)


x11(width=7.5,height=3.3)
# png("Figures/Worldmap_MarineAlienSpecies.png",unit="in",width=7.5,height=3.3,res=300)
layout(matrix(1:2,nc=2),widths=c(0.85,0.15))
op <- par(mar=c(0,0,0,0),las=1,cex=0.9,tck=-0.02,mgp=c(2,0.3,0))

plot(st_geometry(spatial_nspec),col=spatial_nspec$col,lwd=0.5)

plot(1:10,axes=F,type="n")
owncolleg(posx=c(0.1,0.2),posy=c(0.1,0.95),col=cols,cex=1.2,digit=2,zval =c(log10(1),log10(10),log10(50),log10(100),log10(200)) ,zlim=c(0,max(spatial_nspec$nSpec,na.rm=T)))#,">2000"
rect(0.3,-2,1,2,col="white",border="white")
text(x=0.4,y=c(-0.65,0.35,1.05,1.36,1.67),labels=c(1,10,50,100,200))
text(x=0.7,y=0.6,labels="Number of alien taxa",srt=270,cex=1.2)
# text(x=0.8,y=0.6,labels="Number of invasive species",srt=270,cex=1.2)
# mtext(">",side=3,line=-1.27,adj=0.27)
par(op)
# dev.off()



