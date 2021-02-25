##################################################################################
# 
# This script is part of the workflow MaFiR (MArine FIrst Records) to identify
# marine species occurrences in the FirstRecord database.
#
# Point-wise occurrences of species in the FirstRecord database obtained from 
# GBIF and OBIS are matched with terrestrial and marine regions to obtain species
# occurrences per region. 
# 
# Ekin Kaplan, Hanno Seebens, 07.02.2021
##################################################################################


graphics.off()
rm(list=ls())

library(sf)

setwd("/home/hanno/Bioinvasion/InvAccu/MarineExpansion")
# setwd("/scratch/home/hseebens/Bioinvasion/InvAccu/MarineExpansion")


### load data ###############################################################

## First record database with GIBF keys 
firstrecords_GBIF <- read.csv("Data/SpeciesFirstRecords_GBIFkeys.csv",sep=";")
colnames(firstrecords_GBIF)[colnames(firstrecords_GBIF)=="Country"] <- "Region"
colnames(firstrecords_GBIF)[colnames(firstrecords_GBIF)=="NewName"] <- "Species"
uni_spec <- unique(firstrecords_GBIF[,c("Species","speciesKey")])

## (REMOVE CASUALS!!!!)

## (FOR OBIS WE NEED AN ID-SPECIES FILE TO IDENTIFY ALIEN POPULATIONS)
## First record database with OBIS keys
firstrecords_OBIS <- readRDS("Data/RobisCoordinates.rds")

## Polygon file of marine and terrestrial regions
regions <- st_read(dsn="Data/Shapefiles",layer="terraqua",stringsAsFactors = F)

regions$Realm <- NA
regions$Realm[!is.na(regions$ECOREGION)] <- "marine"
regions$Realm[!is.na(regions$Region)] <- "terrestrial"

regions$Region[!is.na(regions$ECOREGION)] <- regions$ECOREGION[!is.na(regions$ECOREGION)]
# regions$Region[!is.na(regions$featurecla)] <- regions$name[!is.na(regions$featurecla)]
regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)

## Terrestrial to marine ecoregion file
neighbours <- read.table("Data/Combined MEOW List_Hanno.csv",sep=";",header=T)
neighbours <- subset(neighbours,Action!="remove")


## Identify region of occurrence for each coordinate entry ######################################

## All taxon-region pairs for the identification of alien populations
TaxonRegionPairs <- paste(firstrecords_GBIF$speciesKey,firstrecords_GBIF$Region,sep="_")

## add marine ecoregions to taxon-region pairs 
marine_firstrecs <- merge(neighbours,firstrecords_GBIF,by="Region")
marine_speckeys <- unique(marine_firstrecs[,c("MEOW","speciesKey")])
meow_firstrecords <- unique(marine_firstrecs[,c("MEOW","Species","speciesKey","FirstRecord","Source")])
TaxonRegionPairs <- c(TaxonRegionPairs,paste(marine_speckeys$speciesKey,marine_speckeys$MEOW,sep="_"))

## list species which are clearly non-marine, clearly marine and clearly freshwater ######
non_marine <- subset(firstrecords_GBIF,Class=="Insecta" 
                     | Phylum=="Tracheophyta"
                     | Phylum=="Anthocerotophyta"
                     | Phylum=="Bryophyta"
                     | Class=="Arachnida"
                     | Class=="Aves"
                     | Class=="Amphibia"
                     | Class=="Mammalia" # not fully correct, but no marine alien mammal known
                     | Habitat_marine=="0")$speciesKey

marine <- subset(firstrecords_GBIF,Habitat_marine=="1")$speciesKey

freshwater <- unique(subset(firstrecords_GBIF,Habitat_freshwater=="1" & Habitat_marine=="0" & Habitat_terrestrial=="0")$speciesKey)
# subset(firstrecords_GBIF,speciesKey%in%freshwater)$Species


## Identify in which region coordinates fall ############################
nchunks <- 21 # set total number of data files, which contain the GBIF records
nsteps <- 10


number_rows <- 10^6

filename <- "/home/hanno/Storage_large/GBIF/FirstRecords_Dec2020/GBIFrecords_Ekinprocessed/GBIFrecords_FRDB_Cleaned_Dec2020.csv"
tot_nrows <- nrow(fread(filename, select = 1L))

steps <- c(seq(0,tot_nrows,by=number_rows),tot_nrows)

for (i in 2:length(steps)){
  
  # Import processed GBIF files
  all_coords <- fread("/home/hanno/Storage_large/GBIF/FirstRecords_Dec2020/GBIFrecords_Ekinprocessed/GBIFrecords_FRDB_Cleaned_Dec2020.csv",
                      # select=c("speciesKey","decimalLatitude","decimalLongitude"),
                      skip=steps[i-1],
                      nrows=number_rows,
                      header=T)
  print(i)  
}


all_counter <- 0
chunk_out <- all_out <- list()
for (i in 1:length(steps)){ # loop over all chunks of GBIF coordinate data

  # Import processed GBIF files
  all_coords <- fread("/home/hanno/Storage_large/GBIF/FirstRecords_Dec2020/GBIFrecords_Ekinprocessed/GBIFrecords_FRDB_Cleaned_Dec2020.csv",
                      # select=c("speciesKey","decimalLatitude","decimalLongitude"),
                      skip=5,
                      nrows=10 ,
                      header=T)
  if (i<22) all_coords <- readRDS(paste0("../../../Storage_large/GBIF/FirstRecords_Dec2020/GBIFrecords_Ekinprocessed/GBIFrecords_Cleaned",i,".rds")) #rdsfile
  if (i==22){
    all_coords <- readRDS("RobisCoordinates.rds") #rdsfile
    colnames(all_coords)[colnames(all_coords)=="id"] <- "speciesKey"
  } 
  
  # Transform to sf object
  coords_sf <- st_as_sf(all_coords,coords=c("decimalLongitude","decimalLatitude"),crs=st_crs(regions))
  
  # Arrange the steps of for loop
  steps <- ceiling(seq(1,nrow(coords_sf),length.out=nsteps))
  
  # Identify region and alien population
  for (j in 1:(length(steps)-1)){# 
    all_counter <- all_counter + 1
    
    print(paste0(round(all_counter/(nsteps*nchunks)*100,2),"%"))
    
    ## identify region of occurrence
    ptspoly <- st_join(coords_sf[steps[j]:steps[j+1],],regions)
    
    ## identify and keep only alien records
    ptspoly$SpeciesRegion <- paste(ptspoly$speciesKey,ptspoly$Region,sep="_")
    ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%TaxonRegionPairs,]
    
    ## export
    coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
    output <- cbind.data.frame(ptspoly_alien$speciesKey,ptspoly_alien$Region,ptspoly_alien$Realm,coords_mat,stringsAsFactors=F) #
    colnames(output) <- c("speciesKey","Region","Realm","Longitude","Latitude")#
    output <- unique(output)
    
    ## remove non-marine species in marine ecoregions and marine species in terrestrial regions
    output <- subset(output,!(output$speciesKey%in%non_marine & output$Realm=="marine"))
    output <- subset(output,!(output$speciesKey%in%marine & output$Realm=="terrestrial"))

    ## remove entries with a very low number of records per region (requires coordinates in the file 'outpout')
    region_records <- as.data.frame(table(output$speciesKey,output$Region),stringsAsFactors = F)
    region_records <- subset(region_records,Freq>0 & Freq<3)
    remove_taxreg <- paste0(region_records$Var1,"_",region_records$Var2)
    output <- subset(output,!(paste0(output$speciesKey,"_",output$Region)%in%remove_taxreg))
    output <- unique(output[,c("speciesKey","Region","Realm")])
    
    ## identify species with the majority of records in terrestrial realm and remove
    realm_spec <- as.matrix(table(output$speciesKey,output$Realm))
    realm_spec_proc <- round(((realm_spec) / rowSums(realm_spec))*100) # percent records per realm
    marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] > 75]
    output <- subset(output,!(speciesKey%in%marinespec & Realm=="terrestrial")) # remove terrestrial records of marine species
    non_marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] <= 75]
    output <- subset(output,!(speciesKey%in%non_marinespec & Realm=="marine")) # remove marine records of non-marine species
    
    # ## test 
    # test_dat <- unique(merge(unique(output[,c("speciesKey","Region","Realm")]),firstrecords_GBIF[,c("speciesKey","Species","Class","Order")],by="speciesKey"))
    # graphics.off()
    # x11(width=12,height=12)
    # plot(st_geometry(regions),xlim=c(0,10),ylim=c(50,60))
    # points(output[which(output$speciesKey==3189846 & output$Region=="North Sea"),c("Longitude","Latitude")],col="black",pch=16)
    # 
    # subset(firstrecords_GBIF,speciesKey==9809222) #>75, no vascular plants, no birds, no insects
    # tab_realm[which(rownames(tab_realm)=="5277297"),]
    
    ## output ###############
    saveRDS(output,paste0("Data/FirstRecords_TerrMarRegions_min3_",i,"_",j,".rds"))
    
    chunk_out[[j]] <- output
  }
  chunk_records <- unique(do.call("rbind",chunk_out))
  
  all_out[[i]] <- chunk_records
  
  ## output ###############
  # saveRDS(chunk_records,paste0("GBIF_First_Records_",i,".rds"))
}
all_records <- do.call("rbind",all_out)
all_records <- unique(all_records)

all_records_spec <- merge(all_records,uni_spec,by="speciesKey",all.x=T)

# subset(all_records_spec,Realm=="marine")

# ## output ###############
saveRDS(all_records_spec,"Data/FirstRecords_TerrMarRegions_min3.rds")
# all_records_spec <- readRDS("Data/FirstRecords_TerrMarRegions_min3.rds")


# x <- 0
# all_out <- list()
# for (i in 1:21){
#   # x <- x + 1
#   # regs_species <- readRDS(paste0("GBIF_First_Records_",i,".rds"))
#   # all_out[[x]] <- regs_species
#   for (j in 1:9){
#     x <- x + 1
#     print(x)
#     regs_species <- readRDS(paste0("FirstRecords_Coords_min5",i,"_",j,".rds"))
#     regs_species <- regs_species[,c("speciesKey","Region","Realm")]
#     regs_species <- regs_species[!duplicated(regs_species),]
#     
#     all_out[[x]] <- regs_species
#   }
# }
# regs_species <- unique(do.call("rbind",all_out))
# # saveRDS(regs_species,"FirstRecords_Coords_min5.rds")
# all_records_spec <- readRDS("FirstRecords_Coords_min5.rds")


## add first records to marine species-regions combination #######################
colnames(all_records_spec)[colnames(all_records_spec)=="Region"] <- "MEOW"
if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
marine_regs_species <- merge(all_records_spec,meow_firstrecords[,-which(colnames(meow_firstrecords)=="Species")],by=c("speciesKey","MEOW"))
marine_regs_species <- subset(marine_regs_species,Realm=="marine")
marine_regspec_fr <- aggregate(FirstRecord ~ MEOW + Species + speciesKey + Realm + Source,data=marine_regs_species,FUN=min)


## add first records to terrestrial species-regions combination #######################
colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Region"
terr_regs_species <- merge(all_records_spec,firstrecords_GBIF[,-which(colnames(firstrecords_GBIF)=="Species")],by=c("speciesKey","Region"))
terr_regs_species <- subset(terr_regs_species,Realm!="marine")
terr_regspec_fr <- aggregate(FirstRecord ~ Region + Species + speciesKey + Realm + Source,data=terr_regs_species,FUN=min)

## combine terrestrial and marine first records #################################
colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Region"
all_regspec_fr <- rbind(marine_regspec_fr,terr_regspec_fr)

## set realm of freshwater species to freshwater ######################
all_regspec_fr$Realm[all_regspec_fr$speciesKey%in%freshwater] <- "freshwater"


## output ###################################################
# write.table(all_regspec_fr,"FirstRecords_TerrMarRegions.csv",sep=";",row.names=F)
write.table(all_regspec_fr,"Data/FirstRecords_TerrMarRegions_min3.csv",sep=";",row.names=F)


# marine_spec <- subset(all_regspec_fr,Region=="North Sea")
# subset(all_regspec_fr,Realm=="freshwater")
# tab_realms <- table(regs_species$Region,regs_species$Realm)
# 
# tab_realms[order(tab_realms[,1]),]
# 
# subset(firstrecords_GBIF,Species=="Sphaerocarpos michelii")

# 
# #Import Alien Records
# fnames <- list.files(pattern = "GBIF_First_Records*")
# xy <- lapply(fnames, readRDS)
# x <- do.call(rbind, xy)
# 
# #Add First Records to Alien Coordinate Records and Save
# lastfile <- merge(x, frmerge, by = c("Taxon", "Region"), all.y = FALSE)
# 
# 
# #Remove coordinates
# lastfile1 <- as.data.frame(cbind(lastfile$Taxon, lastfile$Region, lastfile$FirstRecord))
# lastfile1 <- lastfile1 %>% distinct()
# 
# saveRDS(lastfile1,file = "GBIF_Alien_OCC_FR_1_1.rds")
# 
# #If you want to run a benchmark on how much time is needed to run this code
# time2 <- Sys.time()
# time2 - time1
# 
# #100.000 1.5 minutes, 1.000.000 10 minutes
# 
# 
# 
# 
# 
# #OBIS
# 
# #Clean environment to clear space for OBIS code
# graphics.off()
# rm(list=ls())
# 
# #Set working directory
# setwd("C:/Users/ekink/Desktop/Internship/First Records Project/FirstRecords/Data")
# 
# #Import OBIS occurrence Data
# rdsfile1 <- readRDS("RobisCoordinates")
# 
# #Import Intermediary Files
# setwd("C:/Users/ekink/Desktop/Internship/gbif_downloads/1/one_one")
# sfr_gbif <- read.csv("C:/Users/ekink/Desktop/Internship/gbif_downloads/1/one_one/SpeciesFirstRecords_GBIFkeys.csv", sep=";")
# 
# 
# 
# #Form a Column for Taxons and Regions
# sfr_gbif$TaxonRegion <- paste(sfr_gbif$Species,sfr_gbif$Country,sep="_")
# 
# #Import polygon file
# regions<- st_read(dsn="C:/Users/ekink/Desktop/Internship/Software/GIS/terraqua",layer="terraqua",stringsAsFactors = F)
# 
# #Merge Polygon file with Coordinates
# coords <- st_as_sf(rdsfile1,coords=c("decimalLongitude","decimalLatitude"),
#                    crs=st_crs(regions))
# 
# 
# 
# #Arrange the steps of for loop
# steps <- ceiling(seq(1,nrow(coords),length.out=10))
# 
# 
# counter <- 0
# 
# #Identify and Save alien records
# for (j in 1:(length(steps)-1)){
#         counter <- counter + 1
#         
#         ## identify region of occurrence
#         ptspoly <- st_join(coords[steps[j]:steps[j+1],],regions)
#         
#         ## identify and keep only alien records
#         ptspoly$SpeciesRegion <- paste(ptspoly$scientificName,ptspoly$Region,sep="_")
#         ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%sfr_gbif$TaxonRegion,]
#         
#         ## export
#         coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
#         output <- cbind(ptspoly_alien$scientificName,ptspoly_alien$Region,coords_mat,stringsAsFactors=F)
#         colnames(output) <- c("Taxon","Region","Longitude","Latitude")
#         
#         saveRDS(output,paste0("OBIS_Records_1_",counter,".rds"))}
# 
# #Import Alien Records
# fnames <- list.files(pattern = "OBIS_Records_1_*")
# xy <- lapply(fnames, readRDS)
# x <- do.call(rbind, xy)
# 
# sfr_gbif <- sfr_gbif %>% rename(Taxon = Species, Region = Country)
# sfr_gbif2 <- select(sfr_gbif, Taxon, Region, FirstRecord)
# 
# 
# #Add First Records to Alien Coordinate Records and Save
# lastfile <- merge(x, sfr_gbif2, by = c("Taxon", "Region"), all.y = FALSE)
# 
# #Remove coordinates
# 
# lastfile1 <- as.data.frame(cbind(lastfile$Taxon, lastfile$Region, lastfile$FirstRecord))
# #lastfile1 <- lastfile[!is.na(lastfile$Longitude),]
# #lastfile1 <- lastfile[!is.na(lastfile$Latitude),]
# lastfile1 <- lastfile1 %>% distinct()
# 
# #Save the last file
# saveRDS(lastfile1,file = "OBIS_Alien_OCC_FR_1_1.rds")
# 
# #END
