# run MaFiR workflow

graphics.off()
rm(list=ls())


library(data.table)
library(CoordinateCleaner)
# library(zoo)
library(rgbif)
library(httr)
library(dismo)
library(sf)


### Full path to GBIF download files ###################################
### Note: All files in that folder will be considered as relevant files
path_to_GBIFdownloads <- "/home/hanno/Storage_large/GBIF/SInASdata/GBIF_230221"

name_of_specieslist <- "SInAS_AlienSpeciesDB_2.3.1_FullTaxaList.csv" # has to be stored in Data/Input/ and has to include a column named 'scientificName'

file_name_extension <- "SInAS"



########################################################################
### load scripts #######################################################
source(file.path("R","decompress_file.R")) # a function to decompress large zip files
source(file.path("R","extract_GBIF_columns.R")) # a function to extract from zipped GBIF downloads
source(file.path("R","clean_GBIF_records.R")) # a function to clean GBIF records
# source(file.path("R","request_GBIF_download.R")) # a function to decompress large zip files


########################################################################
### send requests to GBIF ##############################################
n_accounts <- 7
request_GBIF_download(name_of_specieslist)

########################################################################
### extract relevant information from GBIF downloads ###################
extract_GBIF_columns(path_to_GBIFdownloads,file_name_extension)

########################################################################
### clean GBIF records #################################################

## High numbers of records may cause memory issues. The number of records
## can be reduced by setting thin_records to TRUE. Then, coordinates are
## rounded to the second digit and duplicates are removed. For the single
## remaining record at this site, the original (not the rounded) record 
## is kept. Thinning may result in imprecise results when the regions
## considered later in the workflow are small.
thin_records <- T

clean_GBIF_records(path_to_GBIFdownloads,file_name_extension,thin_records)

########################################################################
### identify