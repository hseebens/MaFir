##################################################################################
# Some brief explanation...
#
# This code was written with the intent of requesting the download links of
# large files from GBIF Database to your e-mail inbox.
# 
#
# Ekin Kaplan, Hanno Seebens, 28.12.2020
##################################################################################


graphics.off()
rm(list=ls())
#library(vroom)        #if you need to import very big files
library(rgbif)
library(httr)
library(dismo)
library(splitstackshape)
library(readxl)
library(purrr)
library(dplyr)
library(xlsx)
library(qdap)
library(toxboot)

### Variables #######################################################################

n_chunks <- 21       # number of chunks to divide the species records nearly equally
                     # note that GBIF API only allows 3 simultaneous downloads per account
                     # so, for example, you need to open 7 accounts to simultaneously request
                     # the data in 21 chunks 


### load files #####################################################################

setwd("C:/Users/ekink/Desktop/Internship/First Records Project/FirstRecords/Data")     # Set working directory
#TrialFile <- read_excel("TrialFile.xlsx")
TrialFile <- read.csv2("IntroDat_291120.csv")                           # Load the data table with the species name


names(TrialFile)[names(TrialFile) == "NewName"] <- "Spp"                # Change 'NewName' with your column name that  
                                                                        # contains species names.



TrialFile <- TrialFile %>% distinct(Spp, .keep_all= TRUE)               # Get unique species names


bone <- TrialFile$Spp                                                   # Get the column with your required species names


### get backbone ID's of your species ##############################################


firstlist <- list()
secondlist <- list()
thirdlist <- list()
for (i in 1:length(bone)){
  xtf <- name_backbone(bone[i], limit = 1000, strict = TRUE)            # overrides the max limit
  wthlist <- xtf$speciesKey                                             # max limit= 1000 | get taxonID
  zthlist <- xtf$canonicalName                                          # get also the canonical name
  xthlist <- xtf$matchType                                              # get match type to control if they are "exact"
  firstlist[[i]] <- wthlist
  secondlist[[i]] <- zthlist
  thirdlist[[i]] <- xthlist
if (i%%1000==0) print(i)
}



### Clean your backbone data ###################################################

#thirdlist <- thirdlist[ thirdlist != "NONE" ]
#secondlist <- secondlist[ secondlist != "NULL"]
#firstlist <- firstlist[ firstlist != "NULL"]


firstlist1 <- modify_if(firstlist, is.null, ~c(0))                      # remove null entries
secondlist1 <- modify_if(secondlist, is.null, ~c(0))
thirdlist1 <- modify_if(thirdlist, is.null, ~c(0))


### Merge your bacbone data #####################################################

firstlist1 <- unlist(firstlist1)
secondlist1 <- unlist(secondlist1)
thirdlist1 <- unlist(thirdlist)
firstlist1 <- data.frame(firstlist1)
secondlist1 <- data.frame(secondlist1)
thirdlist1 <- data.frame(thirdlist1)

firstdt <- data.frame(firstlist1$firstlist1, secondlist1$secondlist1, thirdlist1$thirdlist1)
#firstdt <- data.frame(secondlist1$secondlist1, thirdlist1$thirdlist1)

#firstdt <- firstdt[!grepl("HIGHERRANK",firstdt$thirdlist.thirdlist),]
#firstdt <- firstdt[!grepl("FUZZY",firstdt$thirdlist.thirdlist),] 

#secondt <- data.frame(firstlist$firstlist,firstdt$secondlist.secondlist, firstdt$thirdlist.thirdlist)

names(firstdt)[names(firstdt) == "firstlist1.firstlist1"] <- "ID"                   # Change column names
names(firstdt)[names(firstdt) == "secondlist1.secondlist1"] <- "Species Name"
names(firstdt)[names(firstdt) == "thirdlist1.thirdlist1"] <- "Match Type"


### Compile your backbone data ####################################################

firstdt <- data.frame(TrialFile$Spp, firstlist1$firstlist1, secondlist1$secondlist1, thirdlist1$thirdlist1)
write.csv2(secondt, "FRDB.csv")

secondt <- firstdt[!firstdt$ID == 0, ]



secondt <- secondt %>% distinct(`ID`, .keep_all= TRUE)

thirdt <- data.frame(secondt$`Species Name`)

names(thirdt)[names(thirdt) == "secondt..Species.Name."] <- "Spp"

uniquegbif2 <- cSplit(thirdt, 'Spp', sep=" ", type.convert=FALSE)


### Learn the number of occurrences per species ####################################


occurencelist <- list()
for (i in 1:length((uniquegbif2$Spp_1))){
  DaLoop <- try(gbif(uniquegbif2$Spp_1[i], uniquegbif2$Spp_2[i], 
                     args=NULL, geo=TRUE, removeZeros=TRUE, 
                     download=FALSE, start=1, end=Inf))                               #find how many occurrences
  occurencelist[[i]] <- DaLoop}

fourthlist <- unlist(occurencelist)
fourthlist <- data.frame(fourthlist)

fourthdt <- data.frame(secondt$`Species Name`, secondt$ID, secondt$`Match Type`, fourthlist$fourthlist)

names(fourthdt)[names(fourthdt) == "secondt..Species.Name."] <- "Species Name"
names(fourthdt)[names(fourthdt) == "secondt.ID"] <- "ID"
names(fourthdt)[names(fourthdt) == "secondt..Match.Type."] <- "Match Type"
names(fourthdt)[names(fourthdt) == "Values"] <- "Value"


### Divide species names into number of chunks (n_chunks) that you desire,     ###################
### based on their occurrence records                                          ###################

ekoam <- with(fourthdt, split(
  ID, cumsum(c(0, head(purrr::accumulate(Value, ~if ((s <- .x + .y) > sum(fourthdt$Value)/ 21) 0 
                                         else s), -1L)) == 0))) 
counter <- 0
x <- 1


### Prepare the requests for GBIF API #################################################

for (j in ekoam) { 
  
  counter <- counter + 1                                         # counts the loop
  secondlist <- j
  
  ## gbif account details (note that x is part of user name and email address)
  ## in this case, we opened accounts and emails such as: 
  ## (ekinhanno1, ekinhanno1@gmail.com), (ekinhanno2, ekinhanno2@gmail.com) and so on for convenience.
  
  user <- paste0("ekinhanno",x)                                  # your gbif.org username 
  pwd <- "seebenskaplan1234"                                     # your gbif.org password (set the same password for all accounts for convenience)
  email <- paste0("ekinhanno",x,"@outlook.com" )                 # your email which you will recieve the download link
  
  
  if (counter %% 3 == 0){                                        # every time counter can be divided by 3,
    x <- x + 1                                                   # set x + 1 => select new GBIF account below. 
    
  }                                                              # note that GBIF API allows up to 
                                                                 # three simultaneous downloads per account.
  
  ## create JSON script for GBIF API #######################################
  
  sub_keys <- paste(secondlist,collapse=",")
  
  json_script <- paste('
                       {
                       "creator": \"',user,'\",
                        "notificationAddresses": [\"',
                          email,'\"
                         ],
                         "sendNotification": true,
                         "format": "SIMPLE_CSV",
                         "predicate": {
                          "type": "and",
                          "predicates": [{
                            "type": "equals",
                            "key": "HAS_COORDINATE",
                            "value": "true"
                          },{
                            "type": "in",
                            "key": "TAXON_KEY",
                            "values": [',sub_keys,']
                          }]
                         }
                       }
                       ',collapse="",sep="")
  
  
  writeLines(json_script,file("GBIFrequest.json"))
  
  
  ### send the reqest 'GBIFrequest.json' to GBIF ##################
  
  url = "http://api.gbif.org/v1/occurrence/download/request"
  
  request <- POST(url = url, 
                  config = authenticate(user, pwd), 
                  add_headers("Content-Type: application/json"),
                  body = upload_file("GBIFrequest.json"), # path to your local file
                  encode = 'json') 
  
  print(request)

}

#END #####