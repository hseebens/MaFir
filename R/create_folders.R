

create_folders <- function (){
  
  ## create output folder #####
  if (!file.exists("Data")){
    dir.create("Data")
  }
  if (!file.exists(file.path("Data","Output"))){
    dir.create(file.path("Data","Output"))
  }
  if (!file.exists(file.path("Data","Output","Intermediate"))){
    dir.create(file.path("Data","Output","Intermediate"))
  }
  if (!file.exists(file.path("Data","Input"))){
    dir.create(file.path("Data","Input"))
  }
  if (!file.exists(file.path("Data","Input","Shapefiles"))){
    dir.create(file.path("Data","Input","Shapefiles"))
  }
}  