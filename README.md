![Header](Results/Figures/ReadMeHeader.png)
#  Overview
This repository includes the Auto_maxent function which is designed to run [MaxEnt](https://biodiversityinformatics.amnh.org/open_source/maxent/) in a user-friendly and mostly unsupervised way. `Auto_maxent()` is a wrapper of the `dismo::maxent` that allows for easier access to MaxEnt internal parameters and adds multiple tools for data preparation, background sampling, model testing and evaluation, and model selection. In addition to `Auto_maxent()` this repository includes other functions that facilitate the pre-processing, data acquisition and filtering, as well as spatial data harmonization, needed to run any species distribution model. Many of these functions are integrated into Auto_maxent but can be found in the `Functions` folder. These functions can be isolated and offer more options than their Auto_maxent integrated versions. Therefore, in this script, you will find mostly everything you need to run some SDMs.

##  The Auto_maxent function:
`Auto_maxent()`has more than XX parameters that control different aspects of model background sampling, presence point selection, variable selection, and model fit. Most can be run using their default values despite the number of parameters. Other parameters such as XX and XX can have a huge impact on model output and therefore some supervision and user input are recommended. Overall we can group the Auto_maxent parameters into three groups:
### Data preparation and variable selection
- `presence_dat`: Presence points for the focus species in `data.frame` or `sf` format. 
- `predictors`: `SpatRaster` with the environmental variables.
- `coords.p`: if `presence_dat` is a `data.frame`, the name of the columns containing the coordinates of the points (long, lat). By default, this parameter is adjusted to work with the [Gbif](https://www.gbif.org/) default `data.frame` output. 
- `min_obs`: Minumun number of observations to run the analysis. If the number of presence points before or after the cleaning, is lower than `min_obs` the analysis stops.
- `rm.dp`: `[TRUE,FALSE]` should duplicated points be removed from the analysis?
- `name.mod`: Prefix for the mods/output routes. By default, it would take the value of the `species` column from a Gbif data set.
- `sp_range`: An `sf` object containing the range of the species or any other spatial unit that we want to use as a template for our study area. By default this value is `NULL` and the study area is drawn by calculating the minimum complex polygon form by `presence_dat`.
- `crs.r`: coordinate reference system (crs) for the points (presence records and background points) as well as for the environmental variables and study area. To speed up the analysis is recommended to set up this parameter beforehand, making sure that all the spatial information is compatible. See the `resample.rast()` function for details about spatial data harmonization. By default this parameter is set to `EPSG:4326` which returns the WGS-84 crs, this is the default crs of Gbif data and helps to transform this information into spatial data.
- `buff_lim`: the buffer, in map/crs units, to draw around the study area. By default this parameter is set to `0`, therefore no buffer is drawn. 
- `n_bk`: Number of background points to sample from the study area. If set to `AUTO`, the number of background points is defined by the product of the number of cells of the environmental `SpatRaster` and `prop_env`. 
- `prop_env`: Proportion of cells to sample when `n_bk` is set to `AUTO`. The default value is `0.25`, therefore the number of background points is equal to a quarter of the total number of cells in our environmental variables.
- `type_bk`: `[Random,BwData,BwData_inv,EnvBK]` Type of background sampling: Random = points are drawn at random from the study area; BwData = points are sampled from the study area based on a [density kernel of the presence data](put reference here); BwData_inv = points are sampled based on the inverse density kernel of the presence data; EnvBK = points are sampled based on a [density kernel build using the first two Principal Components of a PCA build with the environmental data](ref here). 
- `Test_n`: Proportion of points used for model testing. Default the 30% of presence points are used for model testing and model performance evaluation.
- `world_pol`: An `sf` object to refine the study area (e.g. remove unsuitable areas or marine regions). It helps to keep the background points within specific spatial limits.
- `time_macth`: Should we try to coordinate our spatial data spatiotemporal? Right now this is not implemented (the function would return a message).
- `select_var`: 
![Data Preparation]()
### MaxEnt modelling
- `random_features`:
- `beta.val`:
- `n.m`:
- `Mod.route`:
![modelling]()
### Return options
- `return.all`:
