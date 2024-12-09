rm(list=ls())
gc()
#.rs.restartR()
options(java.parameters = "-Xmx5g") # increase the memory space for jave before loading any package
options("rgdal_show_exportToProj4_warnings"="none") # Silence packages updates warnings

setwd(dirname(rstudioapi::getSourceEditorContext()$path)) # Set the working directory to the directory in which the code is stored
td<-tempdir()
dir.create(td,showWarnings = FALSE)
#
#
###//\/\/\/\/><\/\////\/\/\/\/\/\/\///\\\\\\//\/\/\/\////////////////////////></////##-#
##                        Single species SDM pipeline                               ##-#  
###///\/\/\/><\/\/\////\/\/\//\/\/\/\/\/\///\\\\\\\\\..\.\\\\.\\\\\\\\\\\\\\\><\\\\\##-#
#'
#'
# 0. Load the packages----
list.of.packages<-c("sf","terra","tidyverse","geodata","rJava","viridis")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

lapply(list.of.packages,require,character.only=TRUE)
rm(list.of.packages,new.packages)

# 0.1 Some project parameters----
# CRS and spatial standarization
crs_p <- "ESRI:54012"

# 0.1 Load the needed functions----
functions<-"./Functions" %>% list.files(recursive = FALSE,pattern = ".R$",full.names = TRUE)
lapply(functions,function(x) source(x))

# 1. Get the species occurrence data ----
Spatial_spp(sci_sp = "Macropus rufogriseus") # we are going to make a SDM for the red-necked wallaby
sp_points<- list.files("./Data/sp_points",pattern = "Macropus rufogriseus.csv",full.names = TRUE) %>% read.csv()

# a. Species range----
rX <- st_read("./Data/IUCN_sp_Range" %>% list.files(pattern = ".shp",full.names = TRUE)) %>% st_transform(crs = "EPSG:4326")

# General spatial polygons
terra_border <- geodata::world(resolution=3,level=0,path=td) # takes some time
terra_border <- terra_border %>% st_as_sf() %>% st_make_valid() %>% st_transform(crs = "EPSG:4326") 

# country polygon
country_border <- geodata::gadm(country="aus",path=td,level=0)
country_border <- country_border %>% st_as_sf() %>% st_make_valid() %>% st_transform(crs = "EPSG:4326") 

# Export the spatial objects (reduce the potential problems with the objects saved in the temporal folder)
pol_route <- "./Data/GeneralPols" ; dir.create(pol_route,recursive = TRUE,showWarnings = FALSE)

st_write(terra_border,paste(pol_route,"Continent_pol.shp",sep="/"),append=FALSE)
st_write(country_border,paste(pol_route,"Country_pol.shp",sep="/"),append=FALSE)

T_pol <- st_read(paste(pol_route,"Continent_pol.shp",sep="/"))
co_pol <- st_read(paste(pol_route,"Country_pol.shp",sep="/"))       
        
rm(terra_border,country_border)  
gc()   

# 2. Cleaning species spatial data ----
sp_clean <- Prepare_points(points_sp = sp_points,range_sp = rX)
sp_clean %>% write.csv(paste("./Data/sp_points","Clean_points.csv",sep="/"))

# Transform into an sf object
op <- sp_clean %>% st_as_sf(coords=c("decimalLongitude","decimalLatitude"),crs="EPSG:4326")
rm(sp_clean,sp_points) # remove unnecesary objects

# 3. Prepare environmental data ---- 
# 3.1 Dowload Bioclim 2.0 data ----
# We want this to be relatively fast, so lets take the 10sqkm resolution
bio_clim <- "./Data/Environmental variables/BioClim" ; dir.create(bio_clim,recursive = TRUE,showWarnings = FALSE)
geodata::worldclim_country(country = "AUS",res=0.5,path=bio_clim,var="bio")

# list.files(bio_clim,pattern=".tif",full.names = TRUE,recursive = TRUE) %>% rast()

# 3.2 Prepare some land-use information (in this case we are going to use the land-use from the Land-Use Harmonization project, LH1)----
Lu_route<-"./Data/Environmental variables/Land_use" 
dir.create(Lu_route,recursive = TRUE,showWarnings = FALSE)

Lu_vars<-c("trees", "grassland", "shrubs", "cropland", "built", "bare", "snow", "water", "wetland", "mangroves", "moss")
lapply(Lu_vars,function(x) geodata::landcover(var=x,path=Lu_route))

# 3.3 Load the raster information and harmonize the data ----
env_vars <-  c(bio_clim %>% list.files(pattern = ".tif",recursive=TRUE,full.names = TRUE),
               Lu_route %>% list.files(pattern = ".tif",recursive=TRUE,full.names = TRUE))

resample.rast(y=env_vars,results.r="./Data/Environmental variables/study_area") # adapt and export the data

# 3.3.1 Load the resampled environmental layer
env_vars <- "./Data/Environmental variables/study_area" %>% list.files(pattern=".tif",full.names = TRUE,recursive = TRUE) %>% rast()

# 3.3.a Display the spatial information----
#~~~~~~~~~~~~~~~~
dir.create("junk") # Terra is having a hard time right now, this solve the errors with the mask function
#~~~~~~~~~~~~~~~~
# Some parameters
colors <- colorRampPalette(c("#ffe98aff","#efc02fff","#ad9d52ff","#91ad52ff","#008f1dff")) ; colors<-colors(1000)
results_r <- "Results/Figures/" ; results_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

png(paste(results_r,paste0("sp_general_map",".png"),sep="/"),width=15,height = 15,res=600,units="cm")
  par(bg=NA)
  env_vars[["trees"]] %>% terra::crop(co_pol %>% vect()) %>% terra::mask(co_pol %>% vect()) %>% 
    plot(axes=F,col=colors,plg = list(loc = "bottom", size=c(0.5,1), title = "Tree cover"),
         ext=c(112.921112, 159.109222, -45, -9.142176)) 
  
  rX %>% st_geometry() %>% plot(col="tomato" %>% adjustcolor(alpha.f = 0.35),add=TRUE,alpha=0.5,border="tomato3")
  plot(co_pol %>% st_geometry(),add=TRUE)
  plot(op %>% st_geometry(),add=TRUE,pch=19,cex=0.15,col="firebrick4")
  mtext(side=3,adj=0,"Native distribution of Macropus rufogriseus",font=2,cex=1,xpd=TRUE)

dev.off()

# Detailled map of the species distribution
png(paste(results_r,paste0("sp_detailled_map",".png"),sep="/"),width=15,height = 15,res=600,units="cm")
  par(bg=NA)
    env_vars[["trees"]] %>% terra::crop(co_pol %>% vect()) %>% terra::mask(co_pol %>% vect()) %>% 
      plot(axes=F,col=colors,plg = list(loc = "right", size=c(0.5,1),title = "Tree cover"),
           ext=c(138.896599738, 157, -44, -21.978265543))
    
    rX %>% st_geometry() %>% plot(col="tomato" %>% adjustcolor(alpha.f = 0.35),add=TRUE,alpha=0.5,border="tomato3")
    plot(co_pol %>% st_geometry(),add=TRUE)
    plot(op %>% st_geometry(),add=TRUE,pch=19,cex=0.15,col="firebrick4")
   # mtext(side=3,adj=0,"Native distribution of Macropus rufogriseus",font=2,cex=1,xpd=TRUE)

dev.off()
    
# 4. Run MaxEnt for the species and select the best performing models----
# 4.1 Test different parameters for the models ----
models_r <- "Results/Models" ; models_r %>% dir.create(recursive=TRUE,showWarnings = FALSE)

# Selection of the background sample size:
# The number of background points can have a huge impact on how MaxEnt performs, we wanto to have a sample of points that is
# representative of the environmental data, but we also want to keep the number of samples we take to a minimun in order to 
# optimize the running time of the algorithm
#
  plot(env_vars)
  
# Since we have the range of the species we can reduce a bit more the environmental information in order to speed up a bit the
# bk sampling
  red_env_vars <- env_vars %>% crop(rX %>% vect())
  ind <- !names(red_env_vars) %in% c("snow","moss")

# Test the number of background sampling
  # m_vars <- red_env_vars[[names(red_env_vars)[ind]]] %>% as.data.frame()
  # m_vars <- m_vars[complete.cases(m_vars),]
  # 
  # background.sample <- b_sample(x=m_vars[,c(1:5)], min = 20, max=nrow(op)*3,breaks.r = 100, plot.l = T)
  # n_bk <- background.sample$sample.size.bk

# 4.1.a Random background sampling----
r.maxent <- Auto_maxent(presence_dat=op, 
                        predictors=env_vars, 
                        rm.dp = TRUE, 
                        name.mod = "Wallabi_Press", 
                        type_bk = "Random", #[Random,BwData,BwData_inv,EnvBK]
                        world_pol = co_pol, select_var = "NUMERICAL", 
                        sp_range=rX,
                        random_features = FALSE, beta.val = 7, 
                        n_bk = 10000,
                        # Model Selection
                        mod.select = TRUE, n.mods = 5, use.boyce = 0.5
                        )
# Export the results
# Save the rast.files
route_r <- paste(models_r,"r.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(r.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(r.maxent[-c(8,9)],paste(route_r,"r.maxent.rds",sep="/"))

# 4.1.b Presence Sampling background sampling----
pres.maxent <- Auto_maxent(presence_dat=op, predictors=env_vars, min_obs = 1000,
                           rm.dp = TRUE, 
                           name.mod = "Wallabi_Press", 
                           type_bk = "BwData", #[Random,BwData,BwData_inv,EnvBK]
                           world_pol = co_pol, select_var = "NUMERICAL", 
                           sp_range=rX,
                           random_features = FALSE, beta.val = 7, 
                           n_bk = 10000,
                           # Model Selection
                           mod.select = TRUE, n.mods = 5, use.boyce = 0.5
                            )

# Save the rast.files
route_r <- paste(models_r,"pres.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.maxent[-c(8,9)],paste(route_r,"pres.maxent.rds",sep="/"))

# 4.1.c Presence Sampling background sampling----
pres.inv.maxent <- Auto_maxent(presence_dat=op, predictors=env_vars, min_obs = 1000,
                           rm.dp = TRUE, 
                           name.mod = "Wallabi_PressInv", 
                           type_bk = "BwData_inv", #[Random,BwData,BwData_inv,EnvBK]
                           world_pol = co_pol, select_var = "NUMERICAL", 
                           sp_range=rX,
                           random_features = FALSE, beta.val = 7, 
                           n_bk = 10000,
                           # Model Selection
                           mod.select = TRUE, n.mods = 5, use.boyce = 0.5
                            )

# Save the rast.files
route_r <- paste(models_r,"pres.inv.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(pres.inv.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(pres.inv.maxent[-c(8,9)],paste(route_r,"pres.inv.maxent.rds",sep="/"))

# 4.1.d Presence Sampling background sampling----
env.maxent <- Auto_maxent(presence_dat=op, predictors=env_vars, min_obs = 1000,
                           rm.dp = TRUE, 
                           name.mod = "Wallabi_Env", 
                           type_bk = "EnvBK", #[Random,BwData,BwData_inv,EnvBK]
                           world_pol = co_pol, select_var = FALSE, 
                           sp_range=rX,
                           random_features = FALSE, beta.val = 7, 
                           n_bk = 10000,
                           # Model Selection
                           mod.select = TRUE, n.mods = 5, use.boyce = 0.5
                            )

# Save the rast.files
route_r <- paste(models_r,"env.maxent",sep="/") ; dir.create(route_r,recursive=TRUE)
lapply(env.maxent,function(x) if("SpatRaster" %in% class(x)) x %>% writeRaster(paste(route_r,paste0(ifelse(nlyr(x)>1,"mod.preds","avr.preds"),".tif"),sep="/"),overwrite=TRUE))

saveRDS(env.maxent[-c(8,9)],paste(route_r,"env.maxent.rds",sep="/"))

# 5. Compare the distribution of background points----







# 6. Compare model predictions ----
# Parameters for the plot
col_fun<-colorRampPalette(c("grey35","skyblue","orange","firebrick"))

# Combine the different models
avg.mods <- c(env.maxent$avr.preds,pres.inv.maxent$avr.preds,pres.maxent$avr.preds,r.maxent$avr.preds)
names(avg.mods) <- c("Env","InvPres","Pres","Random") # names

panel(avg.mods,col=col_fun(100))

# 7. Predict Wallabi distribution in the UK ----
# We need to collect the environmental data for the UK
# We want this to be relatively fast, so lets take the 10sqkm resolution
bio_clim <- "./Data/Environmental variables/UK/BioClim" ; dir.create(bio_clim,recursive = TRUE,showWarnings = FALSE)
geodata::worldclim_country(country = "GBR",res=0.5,path=bio_clim,var="bio")

Lu_route<-"./Data/Environmental variables/UK/Land_use" ; dir.create(Lu_route,recursive = TRUE,showWarnings = FALSE)
Lu_vars<-c("trees", "grassland", "shrubs", "cropland", "built", "bare", "snow", "water", "wetland", "mangroves", "moss")
lapply(Lu_vars,function(x) geodata::landcover(var=x,path=Lu_route))

# Harmonize the spatial data----
env_vars <-  "./Data/Environmental variables/UK" %>% list.files(pattern = ".tif",recursive=TRUE,full.names = TRUE)
resample.rast(y=env_vars,results.r = "./Data/Environmental variables/UK")

# download the country boundaries----
country_border <- geodata::gadm(country="gbr",path=td,level=0)
country_border <- country_border %>% st_as_sf() %>% st_make_valid() %>% st_transform(crs = "EPSG:4326") 

# load the new environmental information----
new_dat <- "./Data/Environmental variables/UK/Resample_rast.tif" %>% rast()
new_dat <- new_dat %>% mask(country_border %>% vect())

pres.new.pred <- lapply(pres.inv.maxent$mods,predict,new_dat) %>% rast()

# Run the predictions for the different models----








# End of the script
unlink("junk",recursive = TRUE,force = TRUE)
