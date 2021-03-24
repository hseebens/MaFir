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

  # the loaded file is called 'file_downloads'
  load(file=file.path("Data","Output","GBIF_download_requests.RData")) 
  
  for (i in 1:length(file_downloads)){
    occ_download_get(file_downloads[[i]],overwrite=F,path=path_to_GBIFdownloads)
  }
  
  ## alternative if previous does not work, execute the following and copy-paste output to command line
  # for (i in 1:length(file_downloads)){
  #   cat(paste0("wget https://api.gbif.org/v1/occurrence/download/request/",as.vector(file_downloads[[i]]),".zip \n"))
  # }
}
  