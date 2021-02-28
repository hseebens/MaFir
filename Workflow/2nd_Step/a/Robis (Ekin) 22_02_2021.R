## OBIS ######### #Needs R version 3.6 or higher. Ekin had to upgrade to R version 4.0 due 
# to incompability with mapview package.

## install.packages("githubinstall")
## "githubinstall(iobis/robis")


graphics.off()
rm(list=ls())
library(readxl)
library(robis)
library(dplyr)
library(data.table)

setwd("E:/Robis_Downloads")

### import dataset as "RobisTrial"

RobisTrial <- read.csv("C:/Users/ekink/Desktop/Internship/First Records Project/FirstRecords/Data/IntroDat_291120.csv", sep=";")

### If you want to try it first use "RobisTrial <- head(RobisTrial, 500)"

### Change "NewName" with the column name that contains your desired species names!!!

names(RobisTrial)[names(RobisTrial) == "NewName"] <- "Spp"

### Remove Duplicates
lastry <- as.data.frame(sort(unique(RobisTrial$Spp)))

names(lastry)[names(lastry) == "sort(unique(RobisTrial$Spp))"] <- "Spp"

### Download occurrences one by one (since Robis returns an error for 
### big download files.) with a for loop and save as RDS.


z <- 0
for (k in lastry$Spp) {
        
        z = z + 1

        my_occ <- try(occurrence(scientificname = k, 
                         fields = c("scientificName", "decimalLongitude", "decimalLatitude",
                                    "basisOfRecord", "country", "speciesid", "marine")))
        my_occ$spp <- k

        filename <- paste0("ObisData",z,".rds")
        saveRDS(my_occ, file = filename )
}


#### PROCESS OBIS DOWNLOADS ###################


### Get all the individiual occurrence files and remove errors

filelist <- list.files(pattern = "ObisData*")

allfiles <- list()
for (i in filelist) {
  
        Rfile <- readRDS(i)
        
        if(class(Rfile) != "try-error"){
    
                allfiles[[i]] <- Rfile
   
  }
}

df <- rbindlist(allfiles, use.names = TRUE, fill = TRUE)

df <- na.omit(df, cols = c("speciesid", "decimalLongitude", "decimalLatitude", "scientificName", "basisOfRecord"))

df <- unique(df)

### Save the occurrence records as one big file
saveRDS(df, file = "ObisCompleteDownload.rds")


y1 <- readRDS("ObisCompleteDownload.rds")

y2 <- unique(as.data.frame(cbind(y1$spp, y1$speciesid)))

### Get the species codes
saveRDS(y2, file = "OBIS_Codes_and_Species.rds")

#END