rm(list=ls())
gc()
#.rs.restartR()
options(java.parameters = "-Xmx5g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\/\/\/\/\\\\\-\-\-\-\-\-\--\-\-\-\--\/\/\/\/\/\/\##

#             Prepare Environment; Load the stable libraries and the AutoMaxent working space

#/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/////\/\/\/\/\/\\\\\-\-\-\-\-\-\--\-\-\-\--\\/\/\/\.\.\.\##
#
# 1. Install the needed packages to run the SDM pipeline----
install.packages("remotes") ; library("remotes") ; library("fs")

# 1.a Load the stable version of the packages needed to run the SDM pipeline----
Sinfo <- readRDS("./Session_info/SDM_pipeline/Stable.r" %>% list.files(".rds$",full.names = T))

# Older package versions can conflict with newer or older versions. Therefore, we are going to set up a new package library to host the
old.lib <- .libPaths()[1]
new.lib <- "./Libraries/SDM_pipeline" ; new.lib %>% dir.create(showWarnings = FALSE,recursive = TRUE)  

path_home_r(new.lib) # Set up the direction to the new library as the default (first load all the other packages)

# d Install the packages
Packages <- Sinfo$otherPkgs
Pack.nmes <- Sinfo$otherPkgs %>% names()

# Check if packages are already installed
new.packages <- Pack.nmes[!(Pack.nmes %in% installed.packages(lib.loc = new.lib)[,"Package"])]
if(length(new.packages)>0){
  lapply(Packages[names(Packages) %in% new.packages], function(x) install_version(package = x$Package, version = x$Version, upgrade = "never", lib = new.lib))
  
}else{
  print("AutoMaxent dependencies already installed!")  
}

lapply(installed.packages()[,"Package"],require,character.only=TRUE)
rm(Packages,Pack.nmes)


# 2. Load the Needed Packages to run the AutoMaxent function and its related functions----
# 2.0 Load/install the needed packages ----
list.of.packages<-c("httr","tidyverse","remotes")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 2.1. Connect to the AutoMaxent GitHub repository ----
git_hub <- "https://api.github.com/repos/BioDivHealth/AutoMaxent/git/trees/main?recursive=1"
MaxRepo <- GET(git_hub) # Extract the repo information
MaxRepo

# 2.2 Get the route to the functions ----
file_path <- data.frame(unlist(lapply(content(MaxRepo)$tree, function(x) x$path)))
colnames(file_path) = c('Path')
head(file_path)

# Extract routes
file_path <- file_path %>%
  separate(Path,c('folder','filename'),'/') %>%
  filter(folder == 'Functions') %>%
  filter(str_detect(filename,'.R'))

# 2.3. Configure the routes, download, and export scripts ----
raw_route <- "https://raw.githubusercontent.com/BioDivHealth/AutoMaxent/main" #This is the raw route to the gitHub repository
MyRoute <- paste(getwd(),"Functions/AutoMaxent",sep="/")

dir.create(MyRoute,recursive = T,showWarnings = F)

for(i in 1:nrow(file_path)){
  write_lines(content(GET(paste(raw_route,file_path$folder[i],file_path$filename[i],sep="/"))),
              paste(MyRoute,file_path$filename[i],sep="/"))
}

# 4. Configure the environment to run the AutoMaxent function ----
  # a Get the route to the SessionInfo folder
  file_route <- data.frame(unlist(lapply(content(MaxRepo)$tree, function(x) x$path)))
  MySession <- paste("./Session_info/AutoMaxent_session") ; MySession %>% dir.create(recursive = T,showWarnings = F)
  
  # b Download and export the information    
  writeBin(content(GET(paste(raw_route,"Session_info/SessionInfo.rds",sep="/")),"raw"),paste(MySession,"SessionInfo.rds",sep="/"))
  
  # c Load the information and install the needed packages
  Sinfo <- readRDS(MySession %>% list.files(".rds$",full.names = T))

    # Older package versions can conflict with newer or older versions. Therefore, we are going to set up a new package library to host the
    old.lib <- .libPaths()[1]
    new.lib <- "./Libraries/AutoMaxent" ; new.lib %>% dir.create(showWarnings = FALSE,recursive = TRUE)  
    
    # d Install the packages
  Packages <- Sinfo$otherPkgs
  Pack.nmes <- Sinfo$otherPkgs %>% names()

  # Check if packages are already installed
  new.packages <- Pack.nmes[!(Pack.nmes %in% installed.packages(lib.loc = new.lib)[,"Package"])]
  
  if(length(new.packages)>0){
    lapply(Packages[names(Packages) %in% new.packages], function(x) install_version(package = x$Package, 
                                                                                    version = x$Version, upgrade = "never", 
                                                                                    lib = new.lib, force=FALSE))
    
  }else{
    print("AutoMaxent dependencies already installed!")  
  }
 
#   
# End of the Sript
#  