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


get_GBIF_download <- function(path_to_GBIFdownloads){
  
  occ_download_get(file_download,overwrite=F,path=path_to_GBIFdownloads)
  
}
  