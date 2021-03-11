
add_first_records <- function(
  file_name_extension,
  name_of_TaxonLoc
  ){
  
  
  ### GBIF records ##############################################################
  
  ## load data ###############################

  meow_records <- readRDS(file.path("Data","Output","Intermediate","MarineRecords.rds"))
  
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  if (!"eventDate"%in%colnames(SpecRegionData)){
    cat("\n Column 'eventDate' is missing. Cannot add first records. \n")
  }
  
  
  all_records_spec <- readRDS(file.path("Data","Output",paste0("AlienRegions_GBIF_",file_name_extension,".rds")))

  GBIF_keys <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))
  GBIF_keys <- GBIF_keys[,c("speciesKey","scientificName")]
  
  SpecRegionData_keys <- merge(SpecRegionData,GBIF_keys,by="scientificName",all.x=T)

  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Region"] <- "MEOW"
  if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  marine_regs_species <- merge(all_records_spec,meow_records[,-which(colnames(meow_records)=="scientificName")],by=c("speciesKey","MEOW"),all.x=T)
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regs_species$eventDate[is.na(marine_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  marine_regspec_fr <- aggregate(eventDate ~ MEOW + scientificName + speciesKey + Realm,data=marine_regs_species,FUN=min) # + Source
  marine_regspec_fr$eventDate[marine_regspec_fr$eventDate==2500] <- NA # remove dummy variable
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Location"
  terr_regs_species <- merge(unique(all_records_spec[,-which(colnames(all_records_spec)=="scientificName")]),SpecRegionData_keys,by=c("speciesKey","Location"),all.y=T)
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  terr_regs_species$eventDate[is.na(terr_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  terr_regspec_fr <- aggregate(eventDate ~ Location + scientificName + speciesKey + Realm,data=terr_regs_species,FUN=min)# + Source
  terr_regspec_fr$eventDate[terr_regspec_fr$eventDate==2500] <- NA
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Location"
  all_regspec_fr_GBIF <- rbind(marine_regspec_fr,terr_regspec_fr)
  

  
  ### OBIS records ##############################################################
  
  ## load data ###############################
  
  meow_records <- readRDS(file.path("Data","Output","Intermediate","MarineRecords_OBIS.rds"))
  colnames(meow_records)[colnames(meow_records)=="speciesid"] <- "speciesKey"
  
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  if (!"eventDate"%in%colnames(SpecRegionData)){
    cat("\n Column 'eventDate' is missing. Cannot add first records. \n")
  }
  
  all_records_spec <- readRDS(file.path("Data","Output",paste0("AlienRegions_OBIS_",file_name_extension,".rds")))
  colnames(all_records_spec)[colnames(all_records_spec)=="speciesid"] <- "speciesKey"
  
  OBIS_keys <- readRDS(file.path("Data","Output","OBISTaxa_SInAS_110321.rds"))
  OBIS_keys <- OBIS_keys[,c("speciesid","scientificName")] # scientificName in OBIS is differnt to GBIF!
  colnames(OBIS_keys)[colnames(OBIS_keys)=="speciesid"] <- "speciesKey"
  OBIS_keys$speciesKey <- as.numeric(OBIS_keys$speciesKey)

  SpecRegionData_keys <- merge(SpecRegionData,OBIS_keys,by.x="Taxon",by.y="scientificName")
  
  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Location"] <- "MEOW"
  if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  marine_regs_species <- merge(all_records_spec,meow_records[,-which(colnames(meow_records)=="scientificName")],by=c("speciesKey","MEOW"),all.x=T)
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regs_species$eventDate[is.na(marine_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
  marine_regspec_fr <- aggregate(eventDate ~ MEOW + scientificName + speciesKey + Realm,data=marine_regs_species,FUN=min) # + Source
  marine_regspec_fr$eventDate[marine_regspec_fr$eventDate==2500] <- NA # remove dummy variable
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Location"
  terr_regs_species <- merge(unique(all_records_spec[,-which(colnames(all_records_spec)=="scientificName")]),SpecRegionData_keys,by=c("speciesKey","Location"))#,all.y=T
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  if (nrow(terr_regs_species)>0){
    terr_regs_species$eventDate[is.na(terr_regs_species$eventDate)] <- 2500 ## dummy variable to keep records in aggregate
    terr_regspec_fr <- aggregate(eventDate ~ Location + scientificName + speciesKey + Realm,data=terr_regs_species,FUN=min)# + Source
    terr_regspec_fr$eventDate[terr_regspec_fr$eventDate==2500] <- NA
  }
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Location"
  all_regspec_fr_OBIS <- rbind(marine_regspec_fr,terr_regspec_fr)

  all_regspec_fr <- rbind(all_regspec_fr_GBIF,all_regspec_fr_OBIS)
  all_regspec_fr <- all_regspec_fr[,c("Location","scientificName","Realm","eventDate")]
    
  ## output ###################################################
  write.table(all_regspec_fr,file.path("Data","Output",paste0("AlienRegionsFirstRecords_",file_name_extension,".csv")),sep=";",row.names=F)
  # all_regspec_fr <- read.table(file.path("Data","Output",paste0("AlienRegionsFirstRecords_",file_name_extension,".csv")),sep=";",header=T,stringsAsFactors = F)
  # all_regspec_fr <- readRDS(file.path("Data","FirstRecords_TerrMarRegions_min3.rds"))

  return(all_regspec_fr)
}
