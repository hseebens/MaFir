## Clean GBIF records ###########################################

clean_GBIF_records <- function(path_to_GBIFdownloads,file_name_extension,thin_records){

  # thin_records <- T
  
  ## identify files to import (i.e., all files within all sub-directories ending with .rds and with 'GBIFrecords_NUMBER_NUMBER' in name)
  allfiles <- list.files(path_to_GBIFdownloads)
  GBIF_records_files <- allfiles[grepl("\\.rds",allfiles)]

  dat_all <- list()
  for (i in 1:length(GBIF_records_files)){ # 
    
    cat(paste0("\n ",i,": ",GBIF_records_files[i],"\n"))
    
    # load file
    dat_sub <- readRDS(file.path(path_to_GBIFdownloads,GBIF_records_files[i]))
    
    # remove duplicates
    dat_sub <- unique(dat_sub)
    
    # remove non-numeric values
    nonnumeric <- is.na(as.numeric(dat_sub$speciesKey)) | is.na(as.numeric(dat_sub$decimalLatitude)) | is.na(as.numeric(dat_sub$decimalLongitude))
    if (any(nonnumeric)){
      dat_sub <- dat_sub[!nonnumeric,]
      
      dat_sub$speciesKey <- as.numeric(dat_sub$speciesKey)
      dat_sub$decimalLatitude <- as.numeric(dat_sub$decimalLatitude)
      dat_sub$decimalLongitude <- as.numeric(dat_sub$decimalLongitude)
    }
    
    # remove wrong coordinates
    ind <- (dat_sub$decimalLatitude>90 | dat_sub$decimalLatitude< -90) |  (dat_sub$decimalLongitude>180 | dat_sub$decimalLongitude< -180)
    dat_sub <- dat_sub[!ind,]
    
    # remove empty records
    ind <- is.na(dat_sub$speciesKey) | is.na(dat_sub$decimalLatitude) | is.na(dat_sub$decimalLongitude)
    dat_sub <- dat_sub[!ind,]
    
    # remove inprecise coordinates
    ind <- nchar(dat_sub$decimalLatitude)<5
    dat_sub <- dat_sub[!ind,]
    ind <- nchar(dat_sub$decimalLongitude)<5
    dat_sub <- dat_sub[!ind,]
    
    n_split <- 2*10^5 # number of records per individual chunks (roughly)
    
    if (nrow(dat_sub)>n_split){
      
      cat(paste0("\n Large data set! Split into smaller pieces.\n"))

      tab_rec <- cumsum(table(dat_sub$speciesKey))
      groups <- (ceiling(tab_rec/n_split))
      group_lvl <- unique(groups)
      
      for (j in 1:length(group_lvl)){
        
        cat(paste0("\n ",i,": ",GBIF_records_files[i]," ",j,"/",length(group_lvl),"\n "))
        
        spec_groups <- names(groups)[groups==group_lvl[j]]
        dat_sub_sub <- subset(dat_sub,speciesKey%in%spec_groups)
        
        # thin records by removing duplicated rounded coordinates
        if (thin_records){
          dat_thinned <- list()
          
          cat("\n Record thinning is enabled!  \n")
          
          for (k in 1:length(unique(dat_sub_sub$speciesKey))){
            dat_spec <- subset(dat_sub_sub,speciesKey==unique(dat_sub_sub$speciesKey)[k])
            rounded_lat <- round(dat_spec$decimalLatitude,2) # round coordinates for thinning
            rounded_lon <- round(dat_spec$decimalLongitude,2) # round coordinates for thinning
            ind <- duplicated(cbind(rounded_lon,rounded_lat))
            dat_thinned[[k]] <- dat_spec[!ind,]
          }
          dat_sub_sub <- rbindlist(dat_thinned)
        }

        # clean records
        if (any(table(dat_sub_sub$speciesKey)>(10^5))){ # outlier test of clean_coordinate is memory-consuming for species with many records
          
          cat(paste0("\n No outlier test possible for ",i,"_",j))
          
          dat_cleaned <- clean_coordinates(dat_sub_sub, 
                                           lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                           value ="clean",
                                           tests = c("capitals","centroids", "equal", "gbif", "institutions",  # remove 'seas' test from default
                                                     "zeros")
                                           ) # this outlier methods is more robust compared to the default 'quantile'
          
        } else {
          
          dat_cleaned <- clean_coordinates(dat_sub_sub, 
                                           lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                           value ="clean",
                                           tests = c("capitals","centroids", "equal", "gbif", "institutions", "outliers",  # remove 'seas' test from default
                                                     "zeros"),
                                           outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
          
        }
         
        # intermediate saving of file (just for safety, files can be removed if everything works)
        saveRDS(dat_cleaned,file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".rds")))
      }
      
      # collect data stored to disk
      dat_sub_all <- list()
      for (j in 1:length(group_lvl)){
        
        dat_cleaned <- readRDS(file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".rds")))
        dat_sub_all[[j]] <- dat_cleaned
      }
      dat_cleaned <- rbindlist(dat_sub_all)

      # intermediate saving of file (just for safety, files can be removed if everything works)
      saveRDS(dat_cleaned,file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".rds")))
      
      ## remove intermediate files if previous saving was successful
      if (file.exists(file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".rds")))){
        for (j in 1:length(group_lvl)){
          file.remove(file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,"_",j,".rds")))
        }
      }
    } else { # for smaller data sets
      
      # clean records
      dat_cleaned <- clean_coordinates(dat_sub, 
                                       lon = "decimalLongitude", lat = "decimalLatitude", species = "speciesKey", 
                                       value ="clean",
                                       tests = c("capitals","centroids", "equal", "gbif", "institutions", "outliers",  # remove 'seas' test from default
                                                 "zeros"),
                                       # tests = c( "outliers"),
                                       outliers_method = "mad") # this outlier methods is more robust compared to the default 'quantile'
    }

    # intermediate saving of file (just for safety, files can be removed if everything works)
    saveRDS(dat_cleaned,file.path("Data","Intermediate",paste0("GBIFrecords_Cleaned_",file_name_extension,"_",i,".rds")))
    
    dat_all[[i]] <- dat_cleaned
  }
  
  # output
  dat_all_df <- rbindlist(dat_all)
  
  saveRDS(dat_all_df, file.path("Data","Output",paste0("GBIFrecords_Cleaned_All",file_name_extension,".rds")))

}

