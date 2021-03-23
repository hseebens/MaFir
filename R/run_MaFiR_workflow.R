# run MaFiR workflow

graphics.off()
rm(list=ls())




## Note: keep Source in workflow


library(data.table) # for clean_GBIF_records, request_GBIF_download
library(rgbif) # for clean_GBIF_records, request_GBIF_download
library(worrms)
library(robis)
library(CoordinateCleaner) # for clean_GBIF_records
library(httr)
library(dismo)
library(sf)   # for transform_coords_to_regions



###################################################################################
## load functions #################################################################
source(file.path("R","load_functions.R")) # load all required functions



###################################################################################
### Global variables ##############################################################

### Within this workflow, files will be downloaded and stored in these folders
### Note: All files in that folder will be considered as relevant files
path_to_GBIFdownloads <- "/home/hanno/Storage_large/GBIF/FirstRecords_Mar2021"
path_to_OBISdownloads <- "/home/hanno/Storage_large/OBIS"

## has to be stored in Data/Input/ and has to include a column named 'scientificName'
## for taxon names and 'Location' for region names and 'Taxon' (no authority) for habitat check
# name_of_specieslist <- "SInAS_AlienSpeciesDB_2.3.1_FullTaxaList.csv"
name_of_specieslist <- "IntroDat_22Mar2021.csv"
name_of_taxacolumn <- "TaxonName"
name_of_regioncolumn <- "Region"

## Name of file with the information of alien species and regions
name_of_TaxonLoc <- "IntroDat_22Mar2021.csv"

## name of shapefile providing polygons for the new delineation
name_of_shapefile <- "RegionsTerrMarine"

file_name_extension <- "FirstRecords"


###################################################################################
### Check and create folder structure #############################################
create_folders()


###################################################################################
### Obtaining data ################################################################
### send requests to GBIF #########################################################

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

###################################################################################

## send requests to GBIF 
send_GBIF_request(name_of_specieslist,n_accounts,user=user,pwd=pwd,email=email,colname=name_of_taxacolumn)

## get downloads from GBIF (requires running 'send_GBIF_request' first)
get_GBIF_download(path_to_GBIFdownloads)

### extract relevant information from GBIF downloads ##############################
extract_GBIF_columns(path_to_GBIFdownloads,file_name_extension)


###################################################################################
### get OBIS records ##############################################################

get_OBIS_records(name_of_specieslist, path_to_OBISdownloads,colname=name_of_taxacolumn,file_name_extension)


###################################################################################
### Cleaning data #################################################################

## Thinning of coordinates:
## High numbers of records may cause memory issues. The number of records
## can be reduced by setting thin_records to TRUE. Then, coordinates are
## rounded to the second digit and duplicates are removed. For the single
## remaining record at this site, the original (not the rounded) record
## is kept. Thinning may result in imprecise results when the regions
## considered later in the workflow are small.
thin_records <- T

### clean GBIF records ############################################################

clean_GBIF_records(path_to_GBIFdownloads,file_name_extension,thin_records)

### clean OBIS records ############################################################

clean_OBIS_records(path_to_OBISdownloads,file_name_extension,thin_records)
  

###################################################################################
### get alien regions based on coordintates #######################################

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

coords_to_regions_GBIF(name_of_TaxonLoc,
                       name_of_shapefile,
                       name_of_taxacolumn=name_of_taxacolumn,
                       name_of_regioncolumn=name_of_regioncolumn,
                       realm_extension,
                       file_name_extension)

coords_to_regions_OBIS(name_of_TaxonLoc,
                       name_of_shapefile,
                       name_of_taxacolumn=name_of_taxacolumn,
                       name_of_regioncolumn=name_of_regioncolumn,
                       realm_extension,
                       file_name_extension)


########################################################################
## add first records per region (requires 'eventDate' column) ##########
## and produce final output file containing GBIF and OBIS records ######
dat <- add_first_records(file_name_extension,name_of_TaxonLoc,path_to_GBIFdownloads)




