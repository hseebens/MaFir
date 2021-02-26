rm(list=ls())
graphics.off()


library(rgdal)
library(raster)
library(rgeos)
library(sf)

setwd("/home/hanno/Bioinvasion")

## load shapefiles ##############
# regs_shp <- readOGR(dsn="Data/Regions/",layer="RegionsShapefile_OGR_small",stringsAsFactors = F)
regs_shp <- st_read(dsn="Data/Regions/",layer="RegionsShapefile_200121",stringsAsFactors = F)

# meow_shp <- readOGR(dsn="../DATA/Regions/MEOW",layer="meow_ecos",stringsAsFactors = F)
meow_shp <- st_read(dsn="../DATA/Regions/MEOW",layer="meow_ecos",stringsAsFactors = F)

## identify neighbouring (overlapping marine and terrestrial regions) ####
# inters_poly <- st_intersects(regs_shp,st_buffer(meow_shp,1)) # intersect both shapefiles
inters_poly <- st_intersects(regs_shp,meow_shp) # intersect both shapefiles
inters_df <- as.data.frame(inters_poly)
colnames(inters_df) <- c("Region","MarineEcoregion")

## replace row IDs by names
inters_df$Region <- regs_shp$Region[inters_df$Region] 
inters_df$MarineEcoregion <- meow_shp$ECOREGION[inters_df$MarineEcoregion]
# write.table(inters_df,"Data/Regions/MarineTerrestrialNeighbours_buffer1.csv")
# write.table(inters_df,"Data/Regions/MarineTerrestrialNeighbours.csv")

## generate one single spatial object (feature in sf language) (requires a buffer to close minor gaps)
one_regs_poly <- st_union(st_buffer((regs_shp),dist=0.01))
# st_write(one_regs_poly,dsn="Data/Regions/",layer="One_regs_poly",driver="ESRI Shapefile",delete_layer=T)## output


## remove terrestrial part from marine polgons
diff <- st_difference(meow_shp,one_regs_poly)
# diff2 <- st_difference(diff)

# st_write(diff,dsn="Data/Regions/",layer="MEOW_NoOverlap_combined",driver="ESRI Shapefile",delete_layer=T) ## output



x11(width=12)
plot(st_geometry(diff),col="blue")



