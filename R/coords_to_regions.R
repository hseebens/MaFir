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


coords_to_regions <- function(
        name_of_TaxonLoc,
        name_of_shapefile,
        realm_extension,
        file_name_extension
  ){
  
  ### load data ###############################################################
  
  ## get GBIF species keys
  GBIF_specieskeys <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))
  
  ## Taxon list
  SpecNames <-  read.table(file.path("Data","Input",name_of_specieslist),stringsAsFactors = F,header=T)
  
  if (realm_extension  # check if realms should be identified
      & !file.exists(file.path("Data","Output","Intermediate",paste0("Habitats_",name_of_specieslist)))
      & !all(c("Habitat_marine","Habitat_freshwater","Habitat_terrestrial")%in%colnames(SpecNames))){ # check if required columns exist, if not, download data from WoRMS

      cat("\n 'realm_extension==TRUE' requires information about habitats from WoRMS.")
      cat("\n If not provided in the taxon file, it will be obtained now. \n")
      
      SpecNames <- get_WoRMS_habitats(SpecNames) # get habitats for species in WoRMS
      
      write.table(SpecNames,file.path("Data","Output","Intermediate",paste0("Habitats_",name_of_specieslist)))
  }

  ## load taxon list with habitats if existing
  if (realm_extension & file.exists(file.path("Data","Output","Intermediate",paste0("Habitats_",name_of_specieslist)))){
    SpecNames <- read.table(file.path("Data","Output","Intermediate",paste0("Habitats_",name_of_specieslist)),stringsAsFactors = F)
  }

  SpecNames <- merge(SpecNames,GBIF_specieskeys[,c("scientificName","speciesKey")],by="scientificName")  
  
  ## Taxon x region database 
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  ## standardise location names 
  newLocNames <- standardise_location_names(SpecRegionData$Location,file_name_extension)
  if (nrow(newLocNames)!=nrow(SpecRegionData)){
    stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
  } 
  SpecRegionData$Location <- newLocNames$Location
  # ind <- which(SpecRegionData$Location!=newLocNames$Location)
  
  ## get GBIF species keys
  GBIF_specieskeys <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))

  SpecRegionData_keys <- merge(SpecRegionData,GBIF_specieskeys[,c("scientificName","speciesKey")],by="scientificName")  
  uni_spec <- unique(SpecRegionData_keys[,c("scientificName","speciesKey")])

  # ## (FOR OBIS WE NEED AN ID-SPECIES FILE TO IDENTIFY ALIEN POPULATIONS)
  # ## First record database with OBIS keys
  # firstrecords_OBIS <- readRDS("Data/RobisCoordinates.rds")
  
  ## Polygon file of marine and terrestrial regions
  regions <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  regions$ECOREGION[!is.na(regions$ECOREGION)] <- paste(regions$ECOREGION[!is.na(regions$ECOREGION)],"MEOW",sep="_")
  regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)
  
  if (realm_extension){
    regions$Realm <- NA
    regions$Realm[!is.na(regions$ECOREGION)] <- "marine"
    regions$Realm[!is.na(regions$Region)] <- "terrestrial"
    
    regions$Region[!is.na(regions$ECOREGION)] <- regions$ECOREGION[!is.na(regions$ECOREGION)]
    # regions$Region[!is.na(regions$featurecla)] <- regions$name[!is.na(regions$featurecla)]
    
    ## Terrestrial to marine ecoregion file
    neighbours <- read.table(file.path("Data","Input","Combined MEOW List_Hanno.csv"),sep=";",header=T)
    neighbours <- subset(neighbours,Action!="remove")
    neighbours$MEOW <- paste(neighbours$MEOW,"MEOW",sep="_")
    colnames(neighbours) <- c("Location","MEOW","Action") ## ADJUST shapefile AND REMOVE!!!!!!!!!!!!!!!!!!!!
  }
  
  
  ## Identify region of occurrence for each coordinate entry ######################################
  
  ## All taxon-region pairs for the identification of alien populations
  TaxonRegionPairs <- paste(SpecRegionData_keys$speciesKey,SpecRegionData_keys$Location,sep="_")
  
  if (realm_extension){
    
    ## add marine ecoregions to taxon-region pairs 
    marine_terr_recs <- merge(neighbours,SpecRegionData_keys,by="Location")
    marine_speckeys <- unique(marine_terr_recs[,c("MEOW","speciesKey")])
    meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","speciesKey","eventDate")])#,"Source"
    saveRDS(meow_records,file.path("Data","Output","Intermediate","MarineRecords.rds"))
    
    TaxonRegionPairs <- c(TaxonRegionPairs,paste(marine_speckeys$speciesKey,marine_speckeys$MEOW,sep="_"))
  
    ## list species which are clearly non-marine, clearly marine and clearly freshwater ######
    non_marine <- subset(SpecNames,
                           class=="Insecta" 
                         | phylum=="Tracheophyta"
                         | phylum=="Anthocerotophyta"
                         | phylum=="Bryophyta"
                         | class=="Arachnida"
                         | class=="Aves"
                         | class=="Amphibia"
                         | class=="Mammalia" # not fully correct, but no marine alien mammal known
                         | Habitat_marine=="0")$speciesKey
    
    marine <- subset(SpecNames,Habitat_marine=="1")$speciesKey
    
    freshwater <- unique(subset(SpecNames,Habitat_freshwater=="1" & Habitat_marine=="0" & Habitat_terrestrial=="0")$speciesKey)
  }
  
  
  ## Identify in which region coordinates fall ############################
  
  ## check available files in folder 'Intermediate' or 'Output' ##################
  folder <- file.path("Output","Intermediate")
  available_files <- list.files(file.path("Data","Output","Intermediate"))
  available_files <- available_files[grep("GBIFrecords_Cleaned_",available_files)]
  if (length(available_files)==0){
    folder <- "Output"
    available_files <- list.files(file.path("Data","Output"))
    available_files <- available_files[grep("GBIFrecords_Cleaned_",available_files)]
  }
  if (length(available_files)==0){
    cat("\n No available files with coordinates found! \n")
  }
  
  nchunks <- length(available_files) # set total number of data files, which contain the GBIF records
  nsteps <- 10
  
  all_counter <- 0
  chunk_out <- all_out <- list()
  for (i in 1:nchunks){ # loop over all chunks of GBIF coordinate data
    
    # Import processed GBIF files
    all_coords <- readRDS(file.path("Data",folder,available_files[i])) #rdsfile

    # Transform to sf object
    coords_sf <- st_as_sf(all_coords,coords=c("decimalLongitude","decimalLatitude"),crs=st_crs(regions))
    
    # Arrange the steps of for loop to avoid large memory consumption
    steps <- ceiling(seq(1,nrow(coords_sf),length.out=nsteps))
    
    # Identify region and alien population
    for (j in 1:(length(steps)-1)){# 
      all_counter <- all_counter + 1
      
      print(paste0(round(all_counter/((nsteps-1)*nchunks)*100,2),"%"))
      
      ## identify region of occurrence
      ptspoly <- st_join(coords_sf[steps[j]:steps[j+1],],regions)
      
      ## identify and keep only alien records
      ptspoly$SpeciesRegion <- paste(ptspoly$speciesKey,ptspoly$Region,sep="_")
      ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%TaxonRegionPairs,]
      
      ## export
      coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$speciesKey,ptspoly_alien$Region,ptspoly_alien$Realm,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesKey","Region","Realm","Longitude","Latitude")#
      } else { #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$speciesKey,ptspoly_alien$Region,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesKey","Region","Longitude","Latitude")#
      }
      output <- unique(output)
      
      ## remove non-marine species in marine ecoregions and marine species in terrestrial regions
      if (realm_extension){
        output <- subset(output,!(output$speciesKey%in%non_marine & output$Realm=="marine"))
        output <- subset(output,!(output$speciesKey%in%marine & output$Realm=="terrestrial"))
      }
      
      ## remove entries with a very low number of records per region (requires coordinates in the file 'outpout')
      #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
      region_records <- as.data.frame(table(output$speciesKey,output$Region),stringsAsFactors = F)
      region_records <- subset(region_records,Freq>0 & Freq<3)
      remove_taxreg <- paste0(region_records$Var1,"_",region_records$Var2)
      output <- subset(output,!(paste0(output$speciesKey,"_",output$Region)%in%remove_taxreg))
      output <- unique(output[,c("speciesKey","Region","Realm")])
      
      ## identify species with the majority of records in terrestrial realm and remove
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        realm_spec <- as.matrix(table(output$speciesKey,output$Realm))
        realm_spec_proc <- round(((realm_spec) / rowSums(realm_spec))*100) # percent records per realm
        marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] > 75]
        output <- subset(output,!(speciesKey%in%marinespec & Realm=="terrestrial")) # remove terrestrial records of marine species
        non_marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] <= 75]
        output <- subset(output,!(speciesKey%in%non_marinespec & Realm=="marine")) # remove marine records of non-marine species
      }
      
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
      # saveRDS(output,file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,"_",i,"_",j,".rds")))
      
      chunk_out[[j]] <- output
    }
    chunk_records <- unique(do.call("rbind",chunk_out))
    
    all_out[[i]] <- chunk_records
    
    ## output ###############
    saveRDS(chunk_records,file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,"_",i,".rds")))
  }
  all_records <- do.call("rbind",all_out)
  all_records <- unique(all_records)
  
  all_records_spec <- merge(all_records,uni_spec,by="speciesKey",all.x=T)

  ## set realm of freshwater species to freshwater ######################
  all_records_spec$Realm[all_records_spec$speciesKey%in%freshwater] <- "freshwater"
  
  # ## output ###############
  saveRDS(all_records_spec,file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))
  # all_records_spec <- readRDS(file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))
  
  ## remove intermediate files if previous saving was successful
  if (file.exists(file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))){
    for (i in 1:nchunks){ # loop over all chunks of coordinate data
      for (j in 1:(length(steps)-1)){# 
        file.remove(file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,"_",i,"_",j,".rds")))
      }
    }
  }
  # return(all_records_spec)
}

# all_records_spec <-  readRDS(file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))
# 
# x <- 0
# all_out <- list()
# for (i in 1:21){
#   # x <- x + 1
#   # regs_species <- readRDS(paste0("/home/hanno/Bioinvasion/Others/EkinInternship/GBIF_First_Records_",i,".rds"))
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
# all_records_spec <- readRDS("/home/hanno/Bioinvasion/Others/EkinInternship/FirstRecords_Coords_min5.rds")
# 
# 
# ## add first records to marine species-regions combination #######################
# colnames(all_records_spec)[colnames(all_records_spec)=="Region"] <- "MEOW"
# if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
# marine_regs_species <- merge(all_records_spec,meow_records[,-which(colnames(meow_records)=="scientificName")],by=c("speciesKey","MEOW"))
# marine_regs_species <- subset(marine_regs_species,Realm=="marine")
# marine_regspec_fr <- aggregate(FirstRecord ~ MEOW + Species + speciesKey + Realm + Source,data=marine_regs_species,FUN=min)
# 
# 
# ## add first records to terrestrial species-regions combination #######################
# colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Region"
# terr_regs_species <- merge(all_records_spec,SpecRegionData[,-which(colnames(SpecRegionData)=="Species")],by=c("speciesKey","Region"))
# terr_regs_species <- subset(terr_regs_species,Realm!="marine")
# terr_regspec_fr <- aggregate(FirstRecord ~ Region + Species + speciesKey + Realm + Source,data=terr_regs_species,FUN=min)
# 
# ## combine terrestrial and marine first records #################################
# colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Region"
# all_regspec_fr <- rbind(marine_regspec_fr,terr_regspec_fr)
# 
# ## set realm of freshwater species to freshwater ######################
# all_regspec_fr$Realm[all_regspec_fr$speciesKey%in%freshwater] <- "freshwater"
# 
# 
# ## output ###################################################
# # write.table(all_regspec_fr,"FirstRecords_TerrMarRegions.csv",sep=";",row.names=F)
# write.table(all_regspec_fr,"Data/FirstRecords_TerrMarRegions_min3.csv",sep=";",row.names=F)


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
