##################################################################################
# 
#
# This code was written with the intent of merging and cleaning the downloaded
# files of GBIF and OBIS occurrence records.
# 
#
# Ekin Kaplan, Hanno Seebens, 28.12.2020
##################################################################################


#Merge GBIF and OBIS Download Files
graphics.off()
rm(list=ls())


library(data.table)
library(CoordinateCleaner)

## Set working directory #####
setwd("/home/hanno/Storage_large/GBIF/FirstRecords_Dec2020/gbif_downloads")
#setwd("/scratch/home/hseebens/Storage_large/GBIF/FirstRecords_Dec2020/gbif_downloads")

## identify files to import (i.e., all files within all sub-directories 
## ending with .rds and with 'GBIFrecords_NUMBER_NUMBER' in name)
allfiles <- list.files(recursive=T)
GBIF_records_files <- allfiles[grepl("\\.rds",allfiles)]
GBIF_records_files <- GBIF_records_files[grepl("GBIFrecords_[1-9]_[1-9]"
                                               ,GBIF_records_files)]
# GBIF_records_files <- GBIF_records_files[grepl("GBIFrecords_NoDupl"
#,GBIF_records_files)]


dat_all <- list()
for (i in 1:length(GBIF_records_files)){
  
  # load file
  dat_sub <- readRDS(GBIF_records_files[i])

  # remove empty records
  ind <- is.na(dat_sub$speciesKey)
  # ind <- is.na(dat_sub$decimalLatitude)
  # ind <- is.na(dat_sub$decimalLongitude)
  dat_sub <- dat_sub[!ind,]
  
  # remove duplicates
  dat_sub <- unique(dat_sub)

  # clean records
  dat_cleaned <- clean_coordinates(dat_sub, 
                    lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                    value ="clean",
                    tests = c("capitals","centroids", "equal", "gbif", "institutions", "outliers",  # remove 'seas' test from default
                                 "zeros"))

  # intermediate saving of file (just for safety, files can be removed if everything works)
  saveRDS(dat_sub,paste0("GBIFrecords_Cleaned",i,".rds"))
  
  # collect data
  dat_all[[i]] <- dat_sub
}

# output
dat_all_df <- rbindlist(dat_all)

saveRDS(dat_all_df, "GBIFrecords_FRDB_Cleaned_Dec2020.rds")


#END