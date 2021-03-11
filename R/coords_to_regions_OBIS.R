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


coords_to_regions_OBIS <- function(
        name_of_TaxonLoc,
        name_of_shapefile,
        realm_extension,
        file_name_extension
  ){
  
  ### load data ###############################################################
  
  ## get GBIF species keys
  ## ADJUST file name once the OBIS download is done !!!!!!!!!!!!!!
  OBIS_specieskeys <- readRDS(file.path("Data","Input","OBISTaxa_SInAS_110321.rds"))
  OBIS_specieskeys$speciesid <- as.numeric(OBIS_specieskeys$speciesid)
  
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

  SpecNames <- merge(SpecNames,unique(OBIS_specieskeys[,c("scientificName","speciesid")]),by.x="species",by.y="scientificName")  
  
  ## Taxon x region database 
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  ## standardise location names 
  newLocNames <- standardise_location_names(SpecRegionData$Location,file_name_extension)
  if (nrow(newLocNames)!=nrow(SpecRegionData)){
    stop("\n Standardisation of location names went wrong. Check standardise_location_names.R in coords_to_regions.R \n")
  } 
  SpecRegionData$Location <- newLocNames$Location
  # ind <- which(SpecRegionData$Location!=newLocNames$Location)

  SpecRegionData_keys <- merge(SpecRegionData,OBIS_specieskeys[,c("scientificName","speciesid")],by.x="Taxon",by.y="scientificName")  
  uni_spec <- unique(SpecRegionData_keys[,c("scientificName","speciesid")])

  ## Polygon file of marine and terrestrial regions
  # regions2 <- st_read(dsn=file.path("Data","Input","Shapefiles"),layer=name_of_shapefile,stringsAsFactors = F)
  regions <- st_read(dsn="/home/hanno/Bioinvasion/Data/Regions",layer="MEOW_NoOverlap_combined",stringsAsFactors = F)
  colnames(regions)[2] <- "Ecoregion"
  regions$Ecoregion[!is.na(regions$Ecoregion)] <- paste(regions$Ecoregion[!is.na(regions$Ecoregion)],"MEOW",sep="_")
  # regions <- regions[is.na(regions$featurecla),] # remove lakes !!!! (no alien distinction available yet)
  
  if (realm_extension){
    regions$Realm <- NA
    regions$Realm[!is.na(regions$Ecoregion)] <- "marine"
    regions$Realm[!is.na(regions$Region)] <- "terrestrial"
    
    regions$Location <- regions$Region
    regions$Location[!is.na(regions$Ecoregion)] <- regions$Ecoregion[!is.na(regions$Ecoregion)]
    # regions$Region[!is.na(regions$featurecla)] <- regions$name[!is.na(regions$featurecla)]
    
    ## Terrestrial to marine ecoregion file
    neighbours <- read.table(file.path("Data","Input","Combined MEOW List_Hanno.csv"),sep=";",header=T)
    neighbours <- subset(neighbours,Action!="remove")
    neighbours$MEOW <- paste(neighbours$MEOW,"MEOW",sep="_")
    colnames(neighbours) <- c("Location","MEOW","Action") ## ADJUST shapefile AND REMOVE!!!!!!!!!!!!!!!!!!!!
  }

  ## Identify region of occurrence for each coordinate entry ######################################
  
  ## All taxon-region pairs for the identification of alien populations
  TaxonRegionPairs <- paste(SpecRegionData_keys$speciesid,SpecRegionData_keys$Location,sep="_")
  
  if (realm_extension){
    
    ## add marine ecoregions to taxon-region pairs 
    marine_terr_recs <- merge(neighbours,SpecRegionData_keys,by="Location")
    marine_speckeys <- unique(marine_terr_recs[,c("MEOW","speciesid")])
    meow_records <- unique(marine_terr_recs[,c("MEOW","scientificName","speciesid","eventDate")])#,"Source"
    saveRDS(meow_records,file.path("Data","Output","Intermediate","MarineRecords_OBIS.rds"))
    
    TaxonRegionPairs <- c(TaxonRegionPairs,paste(marine_speckeys$speciesid,marine_speckeys$MEOW,sep="_"))
  
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
                         | Habitat_marine=="0")$speciesid
    
    marine <- subset(SpecNames,Habitat_marine=="1")$speciesid
    
    freshwater <- unique(subset(SpecNames,Habitat_freshwater=="1" & Habitat_marine=="0" & Habitat_terrestrial=="0")$speciesid)
  }
  
  
  ## Identify in which region coordinates fall ############################
  
  ## check available files in folder 'Intermediate' or 'Output' ##################
  folder <- file.path("Output","Intermediate")
  available_files <- list.files(file.path("Data","Output","Intermediate"))
  available_files <- available_files[grep("OBISrecords_Cleaned_All_",available_files)]
  if (length(available_files)==0){
    folder <- "Output"
    available_files <- list.files(file.path("Data","Output"))
    available_files <- available_files[grep("OBISrecords_Cleaned_All_",available_files)]
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
    all_coords$speciesid <- as.numeric(all_coords$speciesid)

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
      ptspoly$SpeciesRegion <- paste(ptspoly$speciesid,ptspoly$Location,sep="_")
      ptspoly_alien <- ptspoly[ptspoly$SpeciesRegion%in%TaxonRegionPairs,]
      
      ## export
      coords_mat <- as.data.frame(st_coordinates(ptspoly_alien),stringsAsFactors = F)
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$speciesid,ptspoly_alien$Location,ptspoly_alien$Realm,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesid","Location","Realm","Longitude","Latitude")#
      } else { #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        output <- cbind.data.frame(ptspoly_alien$speciesid,ptspoly_alien$Location,coords_mat,stringsAsFactors=F) #
        colnames(output) <- c("speciesid","Location","Longitude","Latitude")#
      }
      output <- unique(output)
      
      ## remove non-marine species in marine ecoregions and marine species in terrestrial regions
      if (realm_extension){
        output <- subset(output,!(output$speciesid%in%non_marine & output$Realm=="marine"))
        output <- subset(output,!(output$speciesid%in%marine & output$Realm=="terrestrial"))
      }
      
      ## remove entries with a very low number of records per region (requires coordinates in the file 'outpout')
      #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
      region_records <- as.data.frame(table(output$speciesid,output$Location),stringsAsFactors = F)
      region_records <- subset(region_records,Freq>0 & Freq<3)
      remove_taxreg <- paste0(region_records$Var1,"_",region_records$Var2)
      output <- subset(output,!(paste0(output$speciesid,"_",output$Location)%in%remove_taxreg))
      output <- unique(output[,c("speciesid","Location","Realm")])
      
      ## identify species with the majority of records in terrestrial realm and remove
      if (realm_extension){ #### ADJUST COLUMN NAMES AFTER MAKING THE SHAPEFILE CONSISTENT!!!!
        realm_spec <- as.matrix(table(output$speciesid,output$Realm))
        realm_spec_proc <- round(((realm_spec) / rowSums(realm_spec))*100) # percent records per realm
        marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] > 75]
        output <- subset(output,!(speciesid%in%marinespec & Realm=="terrestrial")) # remove terrestrial records of marine species
        non_marinespec <- rownames(realm_spec_proc)[realm_spec_proc[,which(colnames(realm_spec_proc)=="marine")] <= 75]
        output <- subset(output,!(speciesid%in%non_marinespec & Realm=="marine")) # remove marine records of non-marine species
      }
      
      # ## test 
      # test_dat <- unique(merge(unique(output[,c("speciesid","Location","Realm")]),firstrecords_GBIF[,c("speciesid","Species","Class","Order")],by="speciesid"))
      # graphics.off()
      # x11(width=12,height=12)
      # plot(st_geometry(regions),xlim=c(0,10),ylim=c(50,60))
      # points(output[which(output$speciesid==3189846 & output$Location=="North Sea"),c("Longitude","Latitude")],col="black",pch=16)
      # 
      # subset(firstrecords_GBIF,speciesid==9809222) #>75, no vascular plants, no birds, no insects
      # tab_realm[which(rownames(tab_realm)=="5277297"),]
      
      ## output ###############
      # saveRDS(output,file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,"_",i,"_",j,".rds")))
      
      chunk_out[[j]] <- output
    }
    chunk_records <- unique(do.call("rbind",chunk_out))
    
    all_out[[i]] <- chunk_records
    
    # ## output ###############
    # saveRDS(chunk_records,file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,"_",i,".rds")))
  }
  all_records <- do.call("rbind",all_out)
  all_records <- unique(all_records)
  
  all_records_spec <- merge(all_records,uni_spec,by="speciesid",all.x=T)

  ## set realm of freshwater species to freshwater ######################
  all_records_spec$Realm[all_records_spec$speciesid%in%freshwater] <- "freshwater"
  
  # ## output ###############
  saveRDS(all_records_spec,file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".rds")))
  # all_records_spec_OBIS <- readRDS(file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".rds")))
  
  # ## remove intermediate files if previous saving was successful
  # if (file.exists(file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,".rds")))){
  #   for (i in 1:nchunks){ # loop over all chunks of coordinate data
  #     for (j in 1:(length(steps)-1)){# 
  #       file.remove(file.path("Data","Output","Intermediate",paste0("AlienRegions_OBIS_",file_name_extension,"_",i,"_",j,".rds")))
  #     }
  #   }
  # }
  # return(all_records_spec)
}
