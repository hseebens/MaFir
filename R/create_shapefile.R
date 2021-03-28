##################################################################################
# 
# This script is part of the workflow DASCO to Downscale Alien Species Checklists
# using Occurrence records from GBIF and OBIS.
#
# The DASCO workflow has been published as ..., which has to be cited when used.
#
# Shapefiles providing terrestrial and marine polygons are merged to create one
# shapefile including both.
#
# Authors: Hanno Seebens, Ekin Kaplan, 28.03.2021
##################################################################################



create_shapefile <- function(terrestrial_polygons,marine_polygons){
  
  ## load shapefiles ##########################
  ## terretrial ###########
  regs_shp <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=terrestrial_polygons,stringsAsFactors = F)
  # regs_shp <- st_read(dsn="InvAccu/MarineExpansion/MaFir_Workflow/Data/Input/Shapefiles/",layer="SInAS_Locations",stringsAsFactors = F)
  crs_orig <- st_crs(regs_shp)
  
  regs_shp$Realm <- "terrestrial"

  ## marine ###########
  # meow_shp <- readOGR(dsn="../DATA/Regions/MEOW",layer="meow_ecos",stringsAsFactors = F)
  meow_shp <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=marine_polygons,stringsAsFactors = F)
  
  meow_shp$Realm <- "marine"
  
  ## freshwater (not implemented yet) ###########
  # lake <- st_read(dsn=file.path("Data","Input","Shapefiles"), layer = "ne_10m_lakes", stringsAsFactors = F)
  
  
  ## change projection to Robinson for st_join and st_buffer
  crs <- "+proj=robin +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m" # Robinson projection
  meow_shp <- st_wrap_dateline(meow_shp, options = c("WRAPDATELINE=YES"))
  meow_shp <- st_transform(meow_shp,crs)
  regs_shp <- st_wrap_dateline(regs_shp, options = c("WRAPDATELINE=YES"))
  regs_shp <- st_transform(regs_shp,crs)
  
  ## identify neighbouring (overlapping marine and terrestrial regions) #########
  inters_poly <- st_intersects(regs_shp,meow_shp) # intersect both shapefiles
  inters_df <- as.data.frame(inters_poly)
  colnames(inters_df) <- c("Location","MarineEcoregion")
  
  ## replace row IDs by names
  inters_df$Location <- regs_shp$Location[inters_df$Location] 
  inters_df$MarineEcoregion <- meow_shp$ECOREGION[inters_df$MarineEcoregion]
  # write.table(inters_df,file.path("Data","Regions","NeighbouringPolygons.csv"))
  
  ### remove overlaps and merge ##################################################
  ## generate one single spatial object (feature in sf language) (requires a buffer to close minor gaps)
  # one_regs_poly <- st_union(st_buffer((regs_shp),dist=0.05))
  # # st_write(one_regs_poly,dsn="Data/Regions/",layer="One_regs_poly",driver="ESRI Shapefile",delete_layer=T)## output
  one_regs_poly <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer="One_regs_poly",stringsAsFactors = F)
  one_regs_poly <- st_wrap_dateline(one_regs_poly, options = c("WRAPDATELINE=YES"))
  one_regs_poly <- st_transform(one_regs_poly,crs)
  one_regs_poly <- st_buffer(one_regs_poly,dist=0)
  
  meow_shp <- st_buffer(meow_shp,dist=0)
  
  
  ## remove terrestrial part from marine polgons
  marine <- st_difference(meow_shp,one_regs_poly)
  
  # st_write(marine,dsn=file.path("Data","Input","Shapefiles"),layer="MEOW_NoOverlap_combined",driver="ESRI Shapefile",delete_layer=T) ## output
  # marine <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer="MEOW_NoOverlap_combined",stringsAsFactors=F) 
  # marine <- st_wrap_dateline(marine, options = c("WRAPDATELINE=YES"))
  # marine <- st_transform(marine,crs)
  
  
  ### Output ##########################################
  
  marine <- marine[c("ECOREGION","PROVINCE","REALM","Realm")]
  colnames(marine) <- c("Ecoregion","Province","Realm_MEOW","Realm","geometry")
  marine$Region <- NA
  marine$Ecoregion <- paste0(marine$Ecoregion,"_MEOW")
  marine <- marine[,c("Region","Ecoregion","Province","Realm_MEOW","Realm","geometry")]
  
  terr <- regs_shp[c("Region","Realm")]
  terr$Ecoregion <- NA
  terr$Province <- NA
  terr$Realm_MEOW <- NA
  terr <- terr[,c("Region","Ecoregion","Province","Realm_MEOW","Realm","geometry")]
  
  terr_marine <- rbind(terr,marine)
  terr_marine <- st_transform(terr_marine,crs=crs_orig)
  
  terr_marine$Location <- NA
  terr_marine$Location[!is.na(terr_marine$Region)] <- terr_marine$Region[!is.na(terr_marine$Region)]
  terr_marine$Location[!is.na(terr_marine$Ecoregion)] <- terr_marine$Ecoregion[!is.na(terr_marine$Ecoregion)]
  # terr_marine$Realm <- NA
  # terr_marine$Realm[!is.na(terr_marine$Region)] <- "terrestrial"
  # terr_marine$Realm[!is.na(terr_marine$Ecoregion)] <- "marine"
  terr_marine <- terr_marine[,c("Location","Region","Ecoregion","Province","Realm_MEOW","Realm","geometry")]
  
  st_write(terr_marine,dsn=file.path("Data","Input","Shapefiles"),layer="RegionsTerrMarine",driver="ESRI Shapefile",delete_layer=T) ## output
  
  # x11(width=12)
  # plot(st_geometry(terr_marine),col="blue")

}


