
add_first_records <- function(
  file_name_extension,
  name_of_TaxonLoc
  ){
  
  
  ## load data ###############################
  all_records_spec <- readRDS("Data",file.path("Data","Intermediate",paste0("AlienRegions_",file_name_extension,".rds")))
  
  meow_records <- readRDS(file.path("Data","Intermediate","MarineRecords.rds"))
  
  SpecRegionData <-  read.table(file.path("Data","Input",name_of_TaxonLoc),stringsAsFactors = F,header=T)
  
  
  ## add first records to marine species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="Region"] <- "MEOW"
  if (any(grepl("FirstRecord",colnames(all_records_spec)))) all_records_spec <- all_records_spec[,-grep("FirstRecord",colnames(regs_species))]
  marine_regs_species <- merge(all_records_spec,meow_records[,-which(colnames(meow_records)=="Species")],by=c("speciesKey","MEOW"))
  marine_regs_species <- subset(marine_regs_species,Realm=="marine")
  marine_regspec_fr <- aggregate(FirstRecord ~ MEOW + Species + speciesKey + Realm + Source,data=marine_regs_species,FUN=min)
  
  
  ## add first records to terrestrial species-regions combination #######################
  colnames(all_records_spec)[colnames(all_records_spec)=="MEOW"] <- "Region"
  terr_regs_species <- merge(all_records_spec,SpecRegionData[,-which(colnames(SpecRegionData)=="Species")],by=c("speciesKey","Region"))
  terr_regs_species <- subset(terr_regs_species,Realm!="marine")
  terr_regspec_fr <- aggregate(FirstRecord ~ Region + Species + speciesKey + Realm + Source,data=terr_regs_species,FUN=min)
  
  ## combine terrestrial and marine first records #################################
  colnames(marine_regspec_fr)[colnames(marine_regspec_fr)=="MEOW"] <- "Region"
  all_regspec_fr <- rbind(marine_regspec_fr,terr_regspec_fr)
  
  
  ## output ###################################################
  # write.table(all_regspec_fr,"FirstRecords_TerrMarRegions.csv",sep=";",row.names=F)
  write.table(all_regspec_fr,"Data/FirstRecords_TerrMarRegions_min3.csv",sep=";",row.names=F)
  
}
