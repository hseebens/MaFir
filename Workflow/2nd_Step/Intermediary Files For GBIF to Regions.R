##################################################################################
# 
#
# This code was written with the intent of creating intermediary files to use
# for later stages.
# 
# 
# Ekin Kaplan, Hanno Seebens, 28.12.2020
##################################################################################

### Intermediary Files that will be used in GBIF to Regions Code


### One time use for merging speciesKey with First Record Database ########


## Import FRDB.csv file that was created in GBIF Download Request ####

FRDB <- read.csv("C:/Users/ekink/Desktop/Internship/gbif_downloads/one_one/FRDB.csv", sep=";")

##Import First Record Database ########

IntroDat_291120 <- read.csv("C:/Users/ekink/Desktop/Internship/First Records Project/FirstRecords/Data/IntroDat_291120.csv",
                            sep=";")

names(FRDB)[names(FRDB) == "TrialFile.Spp"] <- "NewName"
names(FRDB)[names(FRDB) == "ID"] <- "speciesKey"

mergedfrdb <- merge.data.frame(IntroDat_291120, FRDB, by = "NewName")

write.csv2(mergedfrdb, "SpeciesFirstRecords_GBIFkeys.csv")


### A file for merging coordinates with First Records ##########

frmerge <- mergedfrdb %>% select(speciesKey, Country, FirstRecord)
names(frmerge)[names(frmerge) == "speciesKey"] <- "Taxon"
names(frmerge)[names(frmerge) == "Country"] <- "Region"
saveRDS(frmerge, file = "First_Record_Merge_File.rds")

### Species names and speciesKey ##################
keys_species <- sfr_gbif %>% select(speciesKey, Species.Name)
names(keys_species)[names(keys_species) == "Species.Name"] <- "speciesName"
write.csv2(keys_species, file = "Species_Name_Keys.csv")

#END