# run MaFiR workflow

graphics.off()
rm(list=ls())




## Note: keep Source in workflow


library(data.table) # for clean_GBIF_records, request_GBIF_download
# library(zoo)
library(rgbif) # for clean_GBIF_records, request_GBIF_download
library(CoordinateCleaner) # for clean_GBIF_records
library(httr)
library(dismo)
library(sf)   # for transform_coords_to_regions
library(worrms)


### Full path to GBIF and OBIS download folders ###################################
### Within this workflow, files will be downloaded and stored in these folders
### Note: All files in that folder will be considered as relevant files
path_to_GBIFdownloads <- "/home/hanno/Storage_large/GBIF/SInASdata/GBIF_030321"
path_to_OBISdownloads <- "/home/hanno/Storage_large/OBIS"

## has to be stored in Data/Input/ and has to include a column named 'scientificName'
## for taxon names and 'Location' for region names and 'Taxon' (no authority) for habitat check
name_of_specieslist <- "SInAS_AlienSpeciesDB_2.3.1_FullTaxaList.csv"

## Name of file with the information of alien species and regions
name_of_TaxonLoc <- "SInAS_AlienSpeciesDB_2.3.1.csv"

## name of shapefile providing polygons for the new delineation
name_of_shapefile <- "RegionsTerrMarine"

file_name_extension <- "SInAS_2"



########################################################################
### load scripts #######################################################
source(file.path("R","decompress_file.R")) # a function to decompress large zip files
source(file.path("R","extract_GBIF_columns.R")) # a function to extract from zipped GBIF downloads
source(file.path("R","send_GBIF_request.R")) # a function to decompress large zip files
source(file.path("R","get_WoRMS_habitats.R")) # get habitat information from WoRMS
source(file.path("R","clean_GBIF_records.R")) # a function to clean GBIF records
source(file.path("R","clean_OBIS_records.R")) # create shapefile of marine and terrestrial polygons
source(file.path("R","coords_to_regions_GBIF.R")) # identify region for each coordinate
source(file.path("R","coords_to_regions_OBIS.R")) # identify region for each coordinate
source(file.path("R","add_first_records.R")) # add first records per species and region (if available)
source(file.path("R","standardise_location_names.R")) # standardise location names (for matching with shapefile)
source(file.path("R","create_shapefile.R")) # create shapefile of marine and terrestrial polygons


### Obtaining data #####################################################
### send requests to GBIF ##############################################

## GBIF account details ############
## Note that multiple accounts are required for n_accounts>1.
## The accounts have to numbered x=1...n_accounts, while x is part of
## user name and email address. For example, user name and email should be:
## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on.

n_accounts <- 7

## login details for first account (x=1) (the '1' in user name and email
## address will be replaced be account number)
user <- "ekinhanno1"                                  # your gbif.org username
pwd <- "seebenskaplan1234"                                     # your gbif.org password (set the same password for all accounts for convenience)
email <- "ekinhanno1@outlook.com"                 # your email which you will recieve the download link

###################################

## send requests to GBIF 
send_GBIF_request(name_of_specieslist,n_accounts,user=user,pwd=pwd,email=email)

## get downloads from GBIF (requires running 'send_GBIF_request' first)
get_GBIF_download(path_to_GBIFdownloads)

### extract relevant information from GBIF downloads ###################
extract_GBIF_columns(path_to_GBIFdownloads,file_name_extension)


### get OBIS records ###################################################

# to be done...


### Cleaning data ######################################################
### clean GBIF records #################################################

## High numbers of records may cause memory issues. The number of records
## can be reduced by setting thin_records to TRUE. Then, coordinates are
## rounded to the second digit and duplicates are removed. For the single
## remaining record at this site, the original (not the rounded) record
## is kept. Thinning may result in imprecise results when the regions
## considered later in the workflow are small.
thin_records <- T

clean_GBIF_records(path_to_GBIFdownloads,file_name_extension,thin_records)

### clean OBIS records #################################################

thin_records <- F
clean_OBIS_records(path_to_OBISdownloads,file_name_extension,thin_records)
  

########################################################################
### get alien regions based on coordintates ############################

## Create shapefile of terrestrial and marine polygons
## loads and combines shapefiles and stores the final shapefile in Data/Input/Shapefiles
terrestrial_polygons <- "RegionsShapefile_200121" # name of terrestrial shapefile
marine_polygons <- "meow_ecos" # name of marine shapefile

create_shapefile(terrestrial_polygons,marine_polygons)


## Assign coordinates to different realms (terrestrial, freshwater, marine)
## depending on geographic location and additional tests
realm_extension <- TRUE 

## assigning country checklists to marine polygons
# checklist_to_marine <- TRUE

# Region shapefile requires a consistent structure for marine and terrestrial polygons !!!!!

coords_to_regions_GBIF(name_of_TaxonLoc,name_of_shapefile,realm_extension,file_name_extension)

coords_to_regions_OBIS(name_of_TaxonLoc,name_of_shapefile,realm_extension,file_name_extension)


########################################################################
## add first records per region (requires 'eventDate' column) ##########
## and produce final output file containing GBIF and OBIS records ######
dat <- add_first_records(file_name_extension,name_of_TaxonLoc)




