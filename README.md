![Header](Results/Figures/ReadMeHeader.png)
#  Overview
This repository includes the Auto_maxent function which is designed to run [MaxEnt](https://biodiversityinformatics.amnh.org/open_source/maxent/) in a user-friendly and mostly unsupervised way. `Auto_maxent()` is a wrapper of the `dismo::maxent` that allows for easier access to MaxEnt internal parameters and adds multiple tools for data preparation, background sampling, model testing and evaluation, and model selection. In addition to `Auto_maxent()` this repository includes other functions that facilitate the pre-processing, data acquisition and filtering, as well as spatial data harmonization, needed to run any species distribution model. Many of these functions are integrated into Auto_maxent but can be found in the `Functions` folder. These functions can be isolated and offer more options than their Auto_maxent integrated versions. Therefore, in this script, you will find mostly everything you need to run some SDMs.

##  The Auto_maxent function:
`Auto_maxent()`has more than XX parameters that control different aspects of model background sampling, presence point selection, variable selection, and model fit. Most can be run using their default values despite the number of parameters. Other parameters such as XX and XX can have a huge impact on model output and therefore some supervision and user input are recommended. Overall we can group the Auto_maxent parameters into three groups:
### Data preparation and variable selection
- `presence_dat`: Presence points for the focus species in `data.frame` or `sf` format. 
- `predictors`: `SpatRaster` with the environmental variables
- `coords.p`: if `presence_dat` is a `data.frame`, the name of the columns containing the coordinates of the points (long, lat). By default, this parameter is adjusted to work with the [Gbif](https://www.gbif.org/) default `data.frame` output. 
- `min_obs`: Minumun number of observations to run the analysis. If the number of presence points before or after the cleaning, is lower than `min_obs` the analysis stops.
- `rm.dp`: `[TRUE,FALSE]` should duplicated points be removed from the analysis?
- `name.mod`: Prefix for the mods/output routes. By default, it would take the value of the `species` column from a Gbif data set.
- `sp_range`: An `sf` object containing the range of the species or any other spatial unit that we want to use as a template for our study area. By default this value is `NULL` and the study area is drawn by calculating the minimum complex polygon form by `presence_dat`.
- `crs.r`: 
- `buff_lim`:
- `n_bk`:
- `prop_env`:
- `type_bk`:
- `Test_n`:
- `world_pol`:
- `time_macth`:
- `select_var`:
