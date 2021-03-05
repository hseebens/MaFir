
add_first_records <- function(
  file_name_extension,
  name_of_TaxonLoc
  ){
  
  
  ## load data ###############################
  all_records_spec <- readRDS(file.path("Data","Output","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))
  
  meow_records <- readRDS(file.path("Data","Output","Intermediate","MarineRecords.rds"))
  
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  if (!"eventDate"%in%colnames(SpecRegionData)){
    cat("\n Column 'eventDate' is missing. Cannot add first records. \n")
  }
  
  GBIF_keys <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))
  GBIF_keys <- GBIF_keys[,c("speciesKey","scientificName")]
  
  SpecRegionData_keys <- merge(SpecRegionData,GBIF_keys,by="scientificName",all.x=T)

    
  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Region"] <- "MEOW"
  if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  marine_regs_species <- merge(all_records_spec,meow_records[,-which(colnames(meow_records)=="scientificName")],by=c("speciesKey","MEOW"),all.x=T)
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regspec_fr <- aggregate(eventDate ~ MEOW + scientificName + speciesKey + Realm,data=marine_regs_species,FUN=min) # + Source
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Location"
  terr_regs_species <- merge(unique(all_records_spec[,-which(colnames(all_records_spec)=="scientificName")]),SpecRegionData_keys,by=c("speciesKey","Location"),all.y=T)
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  terr_regspec_fr <- aggregate(eventDate ~ Location + scientificName + speciesKey + Realm,data=terr_regs_species,FUN=min)# + Source
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Location"
  all_regspec_fr <- rbind(marine_regspec_fr,terr_regspec_fr)
  
  
  ## output ###################################################
  write.table(all_regspec_fr,file.path("Data","Output",paste0("AlienRegionsFirstRecords_",file_name_extension,".csv")),sep=";",row.names=F)
  all_regspec_fr <- read.table(file.path("Data","Output",paste0("AlienRegionsFirstRecords_",file_name_extension,".csv")),sep=";",header=T,stringsAsFactors = F)
}
