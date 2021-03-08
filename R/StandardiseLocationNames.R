

StandardiseLocationNames <- function(dat){

  dat <- SpecRegionData$Location
  
  dat <- cbind.data.frame(dat,dat,1:length(dat),stringsAsFactors=F)  
  colnames(dat) <- c("Location","Location_orig","order")

  ## load location table #################################################
  regions <- read.table(file.path("Data","Input","AllLocations.csv"),sep=";",stringsAsFactors = F,header=T)
  regions$keywords <- gsub("\\(","\\\\(",regions$keywords)
  regions$keywords <- gsub("\\)","\\\\)",regions$keywords)
  regions$keywords <- tolower(regions$keywords) # set all to lower case for matching
  regions$keywords[regions$keywords==""] <- NA
  regions$Location_lower <- tolower(regions$Location) # set all to lower case for matching
  
  ## prepare data set ############################################

  # dat <- read.table(file.path("Output","Intermediate",paste0(inputfiles[i])),header=T,stringsAsFactors = F)
  
  dat_match1 <- unique(dat[,c("Location","Location_orig")]) ## use another dat set for region matching to keep the original names
  # dat_match1$order <- 1:nrow(dat_match1)
  dat_match1$Location <- gsub("\\xa0|\\xc2", " ",dat_match1$Location) # replace weird white space with recognised white space
  dat_match1$Location <- gsub("^\\s+|\\s+$", "",dat_match1$Location) # trim leading and trailing whitespace
  dat_match1$Location <- gsub("  ", " ",dat_match1$Location) # turn two spaces into one
  dat_match1$Location <- gsub(" \\(the\\)", "",dat_match1$Location) # remove " (the)" 
  dat_match1$Location_lower <- tolower(dat_match1$Location) # set all to lower case for matching

  ## step 1: match names of 'dat' with region names of 'regions'
  dat_match1 <- merge(dat_match1,regions,by.x="Location_lower",by.y="Location_lower",all.x=T)
  
  ## step 3: match names by using keywords in 'regions
  ind_keys <- which(!is.na(regions$keywords))
  for (j in 1:length(ind_keys)){ # loop over rows with multiple keywords
    if (any(grepl("; ",regions$keywords[ind_keys[j]]))){ # check if multiple keywords provided
      keywords <- unlist(strsplit(regions$keywords[ind_keys[j]],"; "))
    } else {
      keywords <- regions$keywords[ind_keys[j]]
    }
    for (k in 1:length(keywords)){
      ind_match <- grep(keywords[k],dat_match1$Location_lower) 
      if (length(unique(regions$Location[ind_keys[j]]))>1) cat(paste0("    Warning: ",keywords[k],"match multiple location names. Refine keywords!"))
      dat_match1$Location[ind_match] <- regions$Location[ind_keys[j]]
      dat_match1$countryCode[ind_match]             <- regions$countryCode[ind_keys[j]]
    }
  }
  
  ## final merging of both data sets with standardised region names
  dat_match1 <- dat_match1[order(dat_match1$order),]
  # if (!identical(dat_match1$Taxon_orig,dat$Taxon_orig)) stop("Data sets not sorted equally!")
  dat$Location <- dat_match1$Location
  
  dat_regnames <- merge(dat,regions[,c("locationID","Location")],by.x="Location",by.y="Location",all.x=T)
  
  ## remove duplicated entries ##############
  ind <- which(!duplicated(dat_regnames))
  dat_regnames <- dat_regnames[ind,]
  
  # ## keep only earliest first record
  # if (any(colnames(dat_regnames)=="eventDate")){
  #   oo <- order(dat_regnames$Location,dat_regnames$Taxon,dat_regnames$eventDate) # sort ascending order
  #   dat_regnames <- dat_regnames[oo,]
  #   ind <- which(!duplicated(dat_regnames[,c("Location","Taxon")])) # identify duplicates (the first match is not counted, only the subsequent duplicates)
  #   dat_regnames <- dat_regnames[ind,] # delete duplicates
  # }
  
  ## output ###############################################################################
  
  missing <- dat_regnames$Location_orig[is.na(dat_regnames$locationID)]
  
  # if (length(missing)>0){ # export missing country names
  #   write.table(sort(unique(missing)),file.path("Output","Check",paste0("Missing_Locations_",FileInfo[i,"Dataset_brief_name"],".csv")),row.names = F,col.names=F)
  # }
  
  dat_regnames <- dat_regnames[!is.na(dat_regnames$locationID),]
  # write.table(dat_regnames,file.path("Output","Intermediate",paste0("Step3_StandardLocationNames_",FileInfo[i,"Dataset_brief_name"],".csv")),row.names=F)

  reg_names <- vector()
  for (i in 1:length(inputfiles)){
    dat <- read.table(file.path("Output","Intermediate",paste0("Step3_StandardLocationNames_",FileInfo[i,"Dataset_brief_name"],".csv")),stringsAsFactors = F,header=T)
    reg_names <- rbind(reg_names,cbind(dat[,c("Location","Location_orig")],FileInfo[i,1]))
  }
  reg_names <- reg_names[reg_names$Location!=reg_names$Location_orig,] # export only region names deviating from the original
  reg_names <- unique(reg_names[order(reg_names$Location),])
  colnames(reg_names) <- c("Location","Location_orig","origDB")
  
  write.table(reg_names,file.path("Output","Translated_LocationNames.csv"),row.names=F)
}

