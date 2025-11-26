##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\/\#
## Function to retrieve the taxonomic information and download spatial information from Gbif
##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~}{];/|\[;][;];][;][;]#
#
Spatial_spp <- function(sci_sp, # Scientific name of the species from which we want to gather spatial information
                        p.route=paste(getwd(),"Data/sp_points",sep="/"), # Folder to store the spatial information
                        t.route=paste(getwd(),"Data/Sp_info",sep="/"),
                        range_sp=NULL, # do we have a shapefile with the range of the species?
                        start_date=2015, # The initial date for the spatial query [the function would look from that date onwards till present date]
                        end_date=NULL,
                        IUCN_api=NULL,
                        Gbif=FALSE
){
  
  # 0. Packages and dependencies:
  if(is.null(IUCN_api)) stop(print("You need to provide a valid IUCN_api key"))
  
  options(iucn_redlist_key=IUCN_api)
  
  list.of.packages<-c("tidyr","rredlist","taxize","data.table","stringr",
                      "sp","rgbif","raster","data.table","dplyr","doParallel","parallel")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # 0.1 Create the exit folder 
  p.route %>% dir.create(recursive=TRUE,showWarnings = FALSE)
  
  # a. Check species names and extract synonims information
  y_tax <- retrieve_syns(spp_name=sci_sp,IUCN_api=IUCN_api,Gbif=Gbif)
  y_sp <- y_tax$Spp_syn
  
  t.route %>% dir.create(showWarnings = FALSE,recursive = TRUE)
  y_tax$TaxDat %>% write.csv(paste(t.route,paste0(sci_sp,".csv"),sep="/"))
  
  # b. with the list of synomyns, dowload all the abailable spatial information
  Download_gbif(sp_list=y_sp, # (Character) List of species from which to dowload spatial information
               initial_date=start_date, # (Numeric/year) By default the function will dowload 500 records for each month of the year, from the year specified till present
               end_date=end_date,
               exit_route=p.route, # (Character) Route to store the dowloaded information
               area=range_sp, # (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format. A WKT shape written as either
               gadm_codes=NULL, # (character) The gadm id of the area occurrences are desired from. https://gadm.org/.
               # locality=NULL, # If iterating around different polygons, the name of the regions or polygon
               n_records=150000)
  
  # Return the species taxonomy
  return(y_tax$TaxDat)
  
}

##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\#
## Extract the taxonomic information and list of sinonyms for a given specie
##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~#
#

retrieve_syns<-function(spp_name,   # [Character] The species name from which to collect taxonomic information
                        n_times=5,
                        IUCN_api=NULL,# [Numeric] Number of times the search is repeated until data is found, default value = 1
                        Gbif=Gbif,  # [Logical] Should we check Gbif for a taxonomic matching of the species
                        ITIS=ITIS # the ITIS server is not being maintained at the time, this might cause errors [set to FALSE for the time being] 
)
{
  options(iucn_redlist_key=IUCN_api)
  
  # 0. Load the packages
  list.of.packages<-c("tidyr","rredlist","taxize","data.table","stringr",
                      "rgbif","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # a. Check the species names, if the name is binomial it runs the query to collect the taxonomic information
  # a.1 Removes special characters and fix the capitalization of the name
  spp.x <- stringr::str_trim(spp_name) %>% gsub(pattern = "[[:punct:]]", replacement = " ") %>% stringr::str_trim("both") %>% gsub(pattern = "  ",replacement = "") # Remove white spaces at the start and end of the string
  
  # Resolve capitalization
  CapSp <- function(x) {
    s <- strsplit(x, " ") %>% unlist()
    
    if(length(s)>1){
      paste(paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2)),
            paste(tolower(s[-1]),collapse = " "), sep = " ")
      
    }else{
      paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2))
    }
  }
  
  spp.x <- CapSp(spp.x)
  
  # a.2 Check if the name is related with a species/sub-specie or other taxonomic class (Class, Order, Family) 
  #     by analyzing the number of terms of the character string  
  correct_name <- NULL
  t_11 <- 1
  
  while(is.null(correct_name) && t_11 <= n_times) {
    
    try(correct_name<-gna_verifier(names=spp.x)
        ,silent = TRUE)
    t_11 <- t_11 + 1
  }
  rm(t_11)  
  
  # Summarize the results of the name checking and correction
  if (!is.null(correct_name)){
    if(nrow(correct_name)!=0){
      y.d<-cbind(spp.x,correct_name[,colnames(correct_name)%in%
                                      c("dataSourceTitleShort","acceptedNameScore","matchedCanonicalFull")]) #
      
      names(y.d)[1]<-"or_name"
      
      spp.x<-y.d$matchedCanonicalFull  # Use the corrected name for the rest of taxnomic querys
      
      ifelse(y.d$matchedCanonicalFull==y.d$or_name,
             y.d$Status<-"Correct",
             y.d$Status<-"Incorrect")
      
    }} else {
      y.d<-data.frame(or_name=spp.x,
                      matchedCanonicalFull=NA,
                      Status="Not_found",
                      dataSourceTitleShort=NA,
                      acceptedNameScore=NA)
    }
  
  # b.Use the corrected or original name to look for taxonomic data----
  #
  # b.3. Get the basic data from the ITIS ----
  if(ITIS){
  
  # b.3.1 Get TSN a reference number that we are going to need to gather information from ITIS----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  TSN <- NULL
  t_4 <- 1
  
  while(is.null(TSN) && t_4 <= 5){ # n_times){ # Takes some time to retry the connection
    try(TSN <- get_tsn_(spp.x,searchtype = "scientific"),silent=T)
    t_4 <- t_4 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(TSN %>% is.null() | class(TSN)!= "list") stop("ITIS Search is likely down!! Check https://itis.gov/ for more information")
  
  if (is.null(TSN[[1]])| nrow(TSN[[1]])==0) { 
    
    ITIS_Present<-FALSE
    ITIS_id<-NA
    ITIS_name<-NA
    ITIS_is_valid<-NA
    ITIS_Phylum<-NA
    ITIS_Class<-NA
    ITIS_Order<-NA
    ITIS_Family<-NA
    
    ITIS_accept<-NA
    ITIS_syn<-NA
    ITIS_N_syn<-NA
    
    
  } else {
    
    tsn <- TSN[[1]]
    tsn_n <- tsn[tsn$nameUsage=="valid",] # Select only valid species names
    
    # Select only the species (get rid of the subspecies)
    tsn_n$length <- lapply(tsn_n$scientificName,function(x) strsplit(x,split=" ") %>% unlist() %>% length()) %>% unlist()
    tsn_n <- tsn_n[tsn_n$length == 2,]
    
  if(length(tsn_n$tsn)==0){ # if there is no information about the species in ITIS we set the data.set and return it
      
      ITIS_Present<-FALSE
      ITIS_id<-NA
      ITIS_name<-NA
      ITIS_is_valid<-NA
      ITIS_Phylum<-NA
      ITIS_Class<-NA
      ITIS_Order<-NA
      ITIS_Family<-NA
      
      ITIS_accept<-NA
      ITIS_syn<-NA
      ITIS_N_syn<-NA
      
    }else{
      ITIS_Present<-TRUE
      
      if(!spp.x %in% tsn_n$scientificName){
        tsn_n <- rbind(data.table(tsn=NA,scientificName=spp.x,commonNames=NA, nameUsage=NA, length=NA),tsn_n)
        ITIS_is_valid <-FALSE # The original name is not valid into ITIS but a synonym is
      
      }else{
        ITIS_is_valid <-TRUE  
        
        }
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      Syn <- NULL
      t_5 <- 1
      
      while(is.null(Syn) && t_5 <= n_times){
        
        try(Syn <- lapply(c(tsn_n$scientificName) ,function(x) synonyms(x,db="itis",accepted=T,rows=1)),silent=T)
        
        t_5 <- t_5 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      # Extract the list of Synonyms
      extract_syns <- function(x){
                                  if("acc_name" %in% names(x[[1]])){
                                      y <- c(x[[1]]$syn_name,x[[1]]$acc_name)
                                      p <- c(x[[1]]$syn_tsn,x[[1]]$acc_tsn)
                                    } else {
                                      y <- x[[1]]$syn_name
                                      p <- x[[1]]$syn_tsn
                                    }
                                    
                                    return(list(Syn=y,id=p))
                                  }
      
        Syn_n <- lapply(Syn, function(x) extract_syns(x)$Syn) %>% unlist()
        Syn_tsn <- lapply(Syn, function(x) extract_syns(x)$id) %>% unlist()
      
      # Remove the subspecies
        if(is.null(Syn_n)){
          Syn_n<-NA
          Syn_tsn<-NA
        
        }else{
          index_x <- lapply(Syn_n,function(x) strsplit(x,split=" ") %>% unlist() %>% length()) %>% unlist() ==2
          incex_c <- !duplicated(Syn_n)
          
          Syn_n <- Syn_n[index_x & incex_c]
          Syn_tsn <- Syn_tsn[index_x & incex_c]  
        }
        
      
      # Include the species 
    if (length(Syn_n)==0){
        ITIS_N_syn<-NA
        ITIS_syn<-NA
      
        } else {
      
        ITIS_N_syn <- length(Syn_n) # number of sp sinonims, Subsp excluded
        ITIS_syn <- paste(Syn_n,collapse=";") # combine the names into a single string with all the sinonims
      }
      
      # Get the upstream taxonomic information
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      # For this we need the ITIS accepted name and reference number
      
      q<-NULL
      t_13<-1
      
      while(is.null(q) && t_13 <= n_times){
        try( q <- lapply(c(Syn_tsn,tsn$tsn),itis_acceptname,silent=TRUE) %>% rbindlist(),silent=T)
            t_13 <- t_13 + 1
          }
      
      # Remove the accepted sub-species and get the acepted name and tsn id
      if(is.null(q)){
        ITIS_id <- NA
        ITIS_name <- NA
          
      }else{
        if(!all(is.na(q$acceptedname))){
          q <- q[lapply(q$acceptedname,function(x) strsplit(x,split=" ") %>% unlist() %>% length()) %>% unlist() ==2,]
      
        }
      ITIS_id <- unique(q$acceptedtsn) %>% paste(collapse=";")
      ITIS_name <- unique(q$acceptedname) %>% paste(collapse=";")
      
      }
      
      j<-NULL
      t_13<-1
      
    if(!is.null(q)){
        while(is.null(j) && t_13 <= n_times){
        
        try(j <- ITIS_id[1] %>% as.numeric() %>% itis_hierarchy(what="full"),silent=T)
        t_13 <- t_13 + 1
        }
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    if(is.null(j)){
      ITIS_Phylum <- NA # Change from lower to upper chase
      ITIS_Class <- NA
      ITIS_Order <- NA
      ITIS_Family <- NA
      
    } else {   
      ITIS_Phylum <- toupper(as.character(j[j$rankname=="phylum",4])) # Change from lower to upper chase
      ITIS_Class <- toupper(as.character(j[j$rankname=="class",4]))
      ITIS_Order <- toupper(as.character(j[j$rankname=="order",4]))
      ITIS_Family <- toupper(as.character(j[j$rankname=="family",4]))
        }
    }
  }
  
  # b.3.2 Unify the taxonomic information from ITIS---- 
  ITIS_data<-data.frame(ITIS_Present,ITIS_is_valid,ITIS_id,ITIS_name,
                        ITIS_Phylum,ITIS_N_syn,ITIS_syn,
                        ITIS_Class,ITIS_Order,ITIS_Family)
  }else{
    ITIS_data<-data.frame(ITIS_Present = FALSE,
                          ITIS_is_valid = NA,
                          ITIS_id = NA,
                          ITIS_name = NA,
                          ITIS_Phylum = NA,
                          ITIS_N_syn = NA,
                          ITIS_syn = NA,
                          ITIS_Class = NA,
                          ITIS_Order = NA,
                          ITIS_Family = NA)
    }
  
  # b.4 Get the data from GBIF----
  # Should we retrieve synonym information from GBIF?
  if(Gbif==TRUE){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    key_1 <- NULL
    t_6 <- 1
    
    while(is.null(key_1) && t_6 <= n_times){
      
      try(key_1 <- get_gbifid_(sci=spp.x)[[1]],silent = TRUE) # get the taxon key
      t_6 <- t_6 + 1
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    if (length(key_1)==0){
      
      GBIF_Present<-"No"
      GBIF_id<-NA
      
      GBIF_name<-NA
      
      GBIF_Phylum<-NA
      GBIF_Class<-NA
      GBIF_Order<-NA
      GBIF_Family<-NA
      
      GBIF_syn<-NA
      GBIF_N_syn<-NA
      
    } else{
      
      GBIF_id<-paste(key_1[,colnames(key_1) %in% c("acceptedusagekey")] %>% unique(),collapse="-")
      GBIF_name<-paste(unique(key_1[,colnames(key_1) %in% c("canonicalname")]),collapse="-")
      GBIF_Present <- TRUE
      
      GBIF_Phylum<-toupper(paste(unique(key_1[,colnames(key_1) %in% c("phylum")] %>% unlist() %>% unique()),collapse="-"))
      GBIF_Class<-toupper(paste(unique(key_1[,colnames(key_1) %in% c("class")] %>% unlist() %>% unique()),collapse="-"))
      GBIF_Order<-toupper(paste(unique(key_1[,colnames(key_1) %in% c("order")] %>% unlist() %>% unique()),collapse="-"))
      GBIF_Family<-toupper(paste(unique(key_1[,colnames(key_1) %in% c("family")] %>% unlist() %>% unique()),collapse="-"))
      GBIF_syn<-paste(key_1[,colnames(key_1) %in% c("canonicalname","species")] %>% unlist() %>% unique(),collapse=";")
      
      GBIF_N_syn<-length(GBIF_syn %>% strsplit(split=";") %>% unlist())
    }
    
    GBif_data<-data.frame(GBIF_Present,GBIF_id,GBIF_name,
                          GBIF_N_syn,GBIF_syn,GBIF_Phylum,
                          GBIF_Class,GBIF_Order,GBIF_Family)
  }
  
  
  #   b.1. Get the basic data from the IUCN Red List ----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # Check if the species is present in the IUCN 
  sp_list <- unique(spp.x)
  x1<-NULL
  x2<-NULL
  
  # Combine all the species to look for them in the IUCN
                        if(ITIS){
  if(ITIS_data$ITIS_Present){
    # Using all the available names
    x1 <- lapply(c(ITIS_data$ITIS_name,ITIS_data$ITIS_syn),strsplit,split=";") %>% unlist()
    x1 <- unique(c(spp.x,x1))
    }
  }
     
  if(Gbif){
  if(GBif_data$GBIF_Present){
    x2 <- lapply(c(GBif_data$GBIF_name,GBif_data$GBIF_syn),strsplit,split=";") %>% unlist()
    x2 <- unique(c(spp.x,x2))
    }
  }
  
  sp_list <- c(sp_list,x1,x2) %>% unique()
  sp_long <- lapply(sp_list, function(x) strsplit(x,split=" ") %>% unlist())
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  IUCN_check <- list()
  
  for(i in 1:length(sp_long)){
    g <- NULL
    t_11 <- 1
    
    while( is.null(g) && t_11 <= n_times){
      try(g <- rl_species(genus=sp_long[[i]][1],species=sp_long[[i]][2]),silent=T)
      t_11 <- t_11 + 1
    }
  
    IUCN_check[[i]]<-g
    rm(g)
    # print(i)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if(lapply(IUCN_check,is.null)%>% unlist() %>% all()){
    g <- NULL 
      }else{
    
    g <- IUCN_check[!c(lapply(IUCN_check,is.null)%>% unlist())]    
    }

  if (is.null(g)) {  
    IUCN_id  <- NA
    IUCN_name<- NA
    IUCN_Phylum   <- NA
    IUCN_Class    <- NA
    IUCN_Order    <- NA
    IUCN_Family   <- NA
    
    IUCN_Category <- NA
    IUCN_Present <- "No"
    
    IUCN_N_syn <- NA
    IUCN_syn <- NA
    
  } else {
    
    IUCN_id  <- unique(g[[1]][["taxon"]][["sis_id"]])
    
    # IUCN_name <- unique(g$scientific_name)
    IUCN_name <- unique(g[[1]][["taxon"]][["scientific_name"]])
    
    # if (is.null(g$phylum)){ IUCN_Phylum <- NA } else {IUCN_Phylum <- unique(g$phylum) }
    if (length(g[[1]][["taxon"]][["phylum_name"]])==0){ IUCN_Phylum <- NA } else {IUCN_Phylum <- unique(g[[1]][["taxon"]][["phylum_name"]]) }
    
    # if (is.null(g$class)){ IUCN_Class <- NA }else {IUCN_Class <- unique(g$class) }
    if (length(g[[1]][["taxon"]][["class_name"]])==0){ IUCN_Class <- NA }else {IUCN_Class <- unique(g[[1]][["taxon"]][["class_name"]]) }
    
    # if (is.null(g$order)){ IUCN_Order <- NA }else {IUCN_Order <- unique(g$order) }
    if (length(g[[1]][["taxon"]][["order_name"]])==0){ IUCN_Order <- NA }else {IUCN_Order <- unique(g[[1]][["taxon"]][["order_name"]]) }
    
    # if (is.null(g$family)){ IUCN_Family <- NA } else {IUCN_Family <- unique(g$family)}
    if (length(g[[1]][["taxon"]][["family_name"]])==0){ IUCN_Family <- NA } else {IUCN_Family <- unique(g[[1]][["taxon"]][["family_name"]])}
    
    id_assesment <- g[[1]][["assessments"]][["assessment_id"]]
    
    if(length(id_assesment)!=0){
      a <- rl_assessment(id=id_assesment[length(id_assesment)])
      IUCN_Category <- a[["red_list_category"]][["code"]]
      
    }else{
      IUCN_Category <- NA
    }
    
    IUCN_Present <-"Yes"
    
    if(length(g[[1]][["taxon"]][["synonyms"]])!=0){
      IUCN_syn <- paste(paste(g[[1]][["taxon"]][["synonyms"]]$genus_name,g[[1]][["taxon"]][["synonyms"]]$species_name,sep=" "),collapse=";")
      IUCN_N_syn <- nrow(g[[1]][["taxon"]][["synonyms"]]) # number of sp sinonims, Subsp excluded
      
    } else {
      IUCN_N_syn<-NA
      IUCN_syn<-NA
    }
  }
  
  # b.2 Combine the IUCN information----
  IUCN_data<-data.frame(IUCN_Present,IUCN_id,IUCN_name,
                        IUCN_Category,IUCN_N_syn,IUCN_syn,
                        IUCN_Phylum,IUCN_Class,
                        IUCN_Order,IUCN_Family)
  
  # C. Return the Taxonomic information
  if(Gbif==TRUE){
    Tax_dat<-cbind(Or_name=spp.x,IUCN_data,ITIS_data,GBif_data)
    Spp_syn<-c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name",
                                                       "ITIS_syn","GBIF_syn")]) %>% paste(collapse = ";") %>% 
      strsplit(split = ";") %>% unlist() %>% unique()
    
    Spp_syn<-Spp_syn[!c(Spp_syn %>% grepl(pattern = "NA"))]
    Spp_syn<-Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
  }else{
    Tax_dat <- cbind(Or_name=spp.x,IUCN_data,ITIS_data)
    Spp_syn <- c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name","ITIS_syn","GBIF_syn")]) %>% 
                paste(collapse = ";") %>% strsplit(split = ";") %>% unlist() %>% unique()
    
    Spp_syn <- Spp_syn[!c(Spp_syn %>% grepl(pattern = "NA"))]
    Spp_syn <- Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
    
  }
  # Message 
  print("Taxonomic search done")
}

# End of the script

#
#################################################################################
#       dowload species GBIF point distribution for a given set of species      #-//\/\/\/\/\/\
#################################################################################
#
# This script automatically dowload the Gbif data for a list of species on an specified area
#
Download_gbif<-function(sp_list, # (Character) List of species from which to dowload spatial information
                       initial_date, # (Numeric/year) By default the function will dowload 500 records for each month of the year, from the year specified till present
                       end_date=NULL,
                       # n_cores=2, # (Numeric) Number of cores commited to the processing (only when sp_list>1)
                       exit_route, # (Character) Route to store the dowloaded information
                       area=NULL, # (character) Searches for occurrences inside a polygon in Well Known Text (WKT) format. A WKT shape written as either
                       gadm_codes=NULL, # (character) The gadm id of the area occurrences are desired from. https://gadm.org/.
                       # locality=NULL, # If iterating around different polygons, the name of the regions or polygon
                       n_records=150000 # Maximun number of records to retrieve
){
  
  # Load the library
  library("rgbif")
  
  # Set the dates for the searching
  if(is.null(end_date)){
  present_date<-Sys.time() %>% format("%Y") %>% as.numeric()
  years=c(initial_date:as.numeric(present_date))
    }else{
  years=c(initial_date:as.numeric(end_date))
    }
  
  
  # Create the output folder for the data
  if(!dir.exists(exit_route)){dir.create(exit_route,showWarnings = FALSE)}
  
  if(length(sp_list)<2){
    
    for (k in 1:length(years)){
      print(paste0("year=",k))
      
      for (j in 1:12){
        
        points <- NULL
        t_11 <- 1
        while(is.null(points) && t_11 <= 50) {
          Sys.sleep(1)
          points<-try(occ_search(scientificName =sp_list,
                                 hasCoordinate=TRUE,
                                 year=years[k],
                                 month=j,
                                 hasGeospatialIssue=FALSE,
                                 limit=n_records,
                                 geometry=area,
                                 gadmGid=gadm_codes),
                      silent = FALSE)
          t_11 <- t_11 + 1
        }
        rm(t_11)
        # 
        
        if(is.null(points$data)){
          next
        }else{
          if(!exists("y_points")){
            y_points<-points$data
          }else{
            y_points<-data.table::rbindlist(list(y_points,points$data),fill=TRUE)
          }}
        rm(points)
        print(paste(k,j,sep="-"))
      }
    }
    gc()
    
    if(exists("y_points")){
      
      if(is.null(sp_list)){
        sp_list<-paste("All_records",locality)
      }
      
      write.csv(y_points,paste(exit_route,paste0(sp_list,".csv"),sep="/"),row.names = FALSE)
      rm(y_points)
    }
    
  }else{
    # Create a temporal folder to host the information until its combined
    Temp_tax<-paste(getwd(),"Temp_tax",sep="/")
    dir.create(Temp_tax,showWarnings = FALSE)
    
    for(f in 1:length(sp_list)){
      t1<-Sys.time()
      sp<-sp_list[f]
      
      for (k in 1:length(years)){
        print(paste0("year=",k))
        
        for (j in 1:12){
          
          points <- NULL
          t_11 <- 1
          
          while(is.null(points) && t_11 <= 50) {
            Sys.sleep(1)
            points<-try(occ_search(scientificName =sp,
                                   hasCoordinate=TRUE,
                                   year=years[k],
                                   month=j,
                                   hasGeospatialIssue=FALSE,
                                   limit=n_records,
                                   geometry=area,
                                   gadmGid=gadm_codes),
                        silent = FALSE)
            t_11 <- t_11 + 1
          }
          rm(t_11)
          
          if(is.null(points$data)){
            next
          }else{
            if(!exists("y_points")){
              y_points<-points$data
            }else{
              y_points<-data.table::rbindlist(list(y_points,points$data),fill=TRUE)
            }}
          rm(points)
          print(paste(k,j,sep="-"))
        }
      }
      gc()
      
      if(exists("y_points")){
        write.csv(y_points,paste(Temp_tax,paste0(sp,".csv"),sep="/"),row.names = FALSE)
        rm(y_points)
      }
    }
    
    # Combine the files from the temporal folder, export it to the relevant folder and 
    # erase the temporal folder
    points_spp<-lapply(list.files(Temp_tax,pattern = "csv",full.names = TRUE),fread)
    # points_spp<-points_spp %>% rbindlist(fill=TRUE)
    
    LL<-lapply(points_spp,function(x) apply(x,2,as.character)%>% as.data.frame()) # Data frames contain different column classes whihc cause problems when merging
    points_spp<-LL %>% rbindlist(fill=TRUE)
    
    write.csv(points_spp,paste(exit_route,paste0(sp_list[1],".csv"),sep="/"),row.names = FALSE) # Export the point information under the original name
    unlink(Temp_tax,recursive=TRUE) # Erase the temporal data
    
  }
  return("GBIF data downloaded") 
}


