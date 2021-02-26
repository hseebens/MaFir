
# graphics.off()
rm(list=ls())

library(sf)

setwd("/home/hanno/Bioinvasion/InvAccu/MarineExpansion")



### load data ###############################################################

all_regspec_fr <- read.table("Data/FirstRecords_TerrMarRegions.csv",sep=";",header=T,stringsAsFactors = F)
colnames(all_regspec_fr)[which(colnames(all_regspec_fr)=="NewName")] <- "Species"
# all_regspec_fr <- read.table("Data/FirstRecords_TerrMarRegions_min3.csv",sep=";",header=T,stringsAsFactors = F)
all_regspec_fr <- subset(all_regspec_fr,Realm=="marine")

## Polygon file of marine and terrestrial regions
marine_regions <- st_read(dsn="Data/Shapefiles",layer="MEOW_NoOverlap_combined",stringsAsFactors = F)
mar_regs <- cbind.data.frame(marine_regions$ECOREGION,marine_regions$PROVINCE,stringsAsFactors=F)
colnames(mar_regs) <- c("Ecoregion","Province")

marine_fr <- merge(all_regspec_fr,mar_regs,by.x="Region",by.y="Ecoregion",all.x=T)
# marine_fr <- aggregate(FirstRecord ~ Province + Species, data=marine_fr,min)
marine_fr <- aggregate(FirstRecord ~ Region + Species, data=marine_fr,min)
marine_fr$FirstRecord <- round(marine_fr$FirstRecord/10)*10

# time_series <- as.data.frame(table(marine_fr$FirstRecord,marine_fr$Province),stringsAsFactors = F)
# colnames(time_series) <- c("Year","Province","nSpec")
time_series <- as.data.frame(table(marine_fr$FirstRecord,marine_fr$Region),stringsAsFactors = F)
colnames(time_series) <- c("Year","Region","nSpec")
time_series$Year <- as.integer(time_series$Year)

time_series <- subset(time_series,Year>1500)
time_series <- subset(time_series,Year<2005)
time_series <- subset(time_series,nSpec>0)

uni_prov <- sort(unique(time_series$Province))
sub_dat <- subset(time_series,Region=="North Sea")

x11()
plot(sub_dat$Year,sub_dat$nSpec)
lines(smooth.spline(sub_dat$Year,sub_dat$nSpec))
abline(v=2005,col="gray")
