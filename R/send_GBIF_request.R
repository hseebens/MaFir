##################################################################################
# 
# This script is part of the workflow MaFiR (MArine FIrst Records) to identify
# marine species occurrences in the FirstRecord database.
#
# In this script, a (usually very long) species list is prepared to download
# occurrence records from GBIF. Depending on the number of records per taxon
# downloads are split into n_chunks of species to obtain downloads of similar
# sizes. This requires n_chunks/3 accounts at GBIF.
#
# Ekin Kaplan, Hanno Seebens, 23.02.2021
##################################################################################




send_GBIF_request <- function(name_of_specieslist,n_accounts,user=user,pwd=pwd,email=email){
  
  #######################################################################################
  ### Variables #########################################################################
  
  n_chunks <- n_accounts * 3       # number of chunks to divide the species records nearly equally
  # note that GBIF API only allows 3 simultaneous downloads per account
  # so, for example, you need to open 7 accounts to simultaneously request
  # the data in 21 chunks 
  
  
  #######################################################################################
  ## load species list to be downloaded #################################################
  
  SpecNames <-  read.table(file.path("Data","Input",name_of_specieslist),stringsAsFactors = F,header=T)$scientificName
  
  SpecNames <- unique(SpecNames)                                      # Get unique species names
  
  
  
  #######################################################################################
  ### get GBIF keys for all species #####################################################
  
  cat("\n Get GBIF keys for taxa \n")
  
  GBIF_speclist <- list()
  x <- 0
  for (i in 1:length(SpecNames)){# loop over all species
    
    specname <- name_backbone(SpecNames[i], limit = 10, strict = TRUE)      # overrides the max limit to increase speed
    if (all(colnames(specname)!="species")) next
    
    x <- x + 1
    GBIF_speclist[[x]] <- c(specname$speciesKey,specname$scientificName,specname$matchType,SpecNames[i])
    
    if (x%%1000==0) print(x)
  }
  GBIF_species <- as.data.frame(do.call("rbind",GBIF_speclist),stringsAsFactors = F)
  colnames(GBIF_species) <- c("speciesKey","scientificName","matchType","Orig_name")
  
  
  ## save intermediate output ######
  write.csv2(GBIF_species, file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))
  # GBIF_species <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFkeys.csv"))
  
  
  
  
  #######################################################################################
  ### Get the number of GBIF records per species ########################################
  ### to split the chunks into roughly equal pieces for download ########################
  
  cat("\n Get number of records per taxon from GBIF \n")
  
  
  spec_genus_split <- strsplit(GBIF_species$scientificName, " ")
  
  GBIF_species$nRecords <- 0
  for (i in 1:length(spec_genus_split)){
    
    # extract number of occurrences  
    options(warn=-1) # suppress warnings
    nRecords <- try(gbif(spec_genus_split[[i]][1], spec_genus_split[[i]][2], 
                         geo=TRUE, removeZeros=TRUE, 
                         download=FALSE, start=1, end=Inf),silent=T)
    options(warn=0) # set default
    
    if (class(nRecords)=="try-error") next
    
    GBIF_species$nRecords[i] <- nRecords
    
    if (i%%1000==0) print(i)
    
  }
  
  ## save intermediate output #############################
  write.csv2(GBIF_species, file.path("Data","Output","Intermediate","SpeciesGBIFnRecords.csv"))
  # GBIF_species <- read.csv2(file.path("Data","Output","Intermediate","SpeciesGBIFnRecords.csv"))
  
  
  
  #######################################################################################
  ### Split taxon list into n_chunks of similar size (nRecords) for download ############
  ### (i.e., split sum of records of all species into n_chunks)
  
  
  GBIF_species$cumsum <- 0
  GBIF_species$cumsum[1] <- GBIF_species$nRecords[1]
  GBIF_species$group <- 1
  x <- 1
  for (i in 2:nrow(GBIF_species)){ # loop over all species
    if (GBIF_species$cumsum[i-1] > sum(GBIF_species$nRecords)/ 21){ # if $cumsum is larger x, start a new cumulative sum
      x <- x + 1 # counter for the number of groups of species belonging to one chunk
      GBIF_species$cumsum[i] <- GBIF_species$nRecords[i]
    } else { # if not, continue with the former cumulative sum
      GBIF_species$cumsum[i] <- GBIF_species$cumsum[i-1] + GBIF_species$nRecords[i]
    }
    GBIF_species$group[i] <- x # set the group number
  }
  
  #######################################################################################
  ### Prepare the requests for GBIF API #################################################
  ### Two approaches are provided below using either several GBIF accounts or a single
  ### account. The first is less convenient but stable, while the latter is a beta version
  
  ### using various GBIF accounts #######################################################
  user_base <- user
  email_base <- email
  
  counter <- 0
  x <- 1
  file_downloads <- list()
  for (j in unique(GBIF_species$group)) {
    
    counter <- counter + 1                                         # counts the loop
    
    ## gbif account details (note that x is part of user name and email address)
    ## in this case, we opened accounts and emails such as:
    ## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on for convenience.
    
    user <- gsub("1",x,user_base)                  # your gbif.org username
    email <- gsub("1",x,email_base)                # your email which you will recieve the download link

    if (counter %% 3 == 0){                                        # every time counter can be divided by 3,
      x <- x + 1                                                   # set x + 1 => select new GBIF account below.
    }                                                              # note that GBIF API allows up to
    
    ## send query of each chunk to GBIF #######################################
    
    sub_keys <- subset(GBIF_species,group==j)$speciesKey
    
    ## prepare requests for GBIF download (no execution!)
    file_downloads[[j]] <- occ_download(
      pred_in("taxonKey", sub_keys),
      pred("hasCoordinate", TRUE),
      pred("hasGeospatialIssue", FALSE),
      format = "SIMPLE_CSV",
      user=user,pwd=pwd,email=email
    )
  }
  
  save(file_downloads,file=file.path("Data","Output","GBIF_download_requests.RData"))
  
  # ### using one GBIF account and different queries #######################################################
  # ### the rgbif functions are beta versions and may not work as expected! ################################
  # 
  # ## GBIF account details ##############################################################################
  # 
  # user <- paste0("ekinhanno",1)                                  # your gbif.org username
  # pwd <- "seebenskaplan1234"                                     # your gbif.org password (set the same password for all accounts for convenience)
  # email <- paste0("ekinhanno",1,"@outlook.com" )                 # your email which you will recieve the download link
  # 
  # queries <- list()
  # for (j in unique(GBIF_species$group)) {
  # 
  #   ## send query of each chunk to GBIF #######################################
  #   
  #   sub_keys <- subset(GBIF_species,group==j)$speciesKey
  #   
  #   ## prepare requests for GBIF download (no execution!)
  #   queries[[j]] <- occ_download_prep(
  #     pred_in("taxonKey", sub_keys),
  #     pred("hasCoordinate", TRUE),
  #     pred("hasGeospatialIssue", FALSE),
  #     format = "SIMPLE_CSV",
  #     user=user,pwd=pwd,email=email
  #   )
  # }
  # 
  # ## execute requests in sequence
  # out <- occ_download_queue(.list = queries, status_ping = 60)
  # out
  
}

# dat <- readRDS("/home/hanno/Storage_large/GBIF/SInASdata/GBIFrecords_Cleaned_AllSInAS.rds")
