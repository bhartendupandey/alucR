alucR - allocation of land use change 
---
# alucR
alucR - Project is a first step to implement a Land Use/Cover Change Model in R (http://www.r-project.org). 
We have been following the basic framework provided by Verburg et al. (2002). The code uses basic R-language and packages and is fully documented. This makes
it possible to easily adapt the code to the users specific needs. 
The large community of R-users might additional help to implement improvements in performance and flexibility to the current procedures.

## Function definition:
aluc(lc, suit, spatial, demand, elas, traj, nochange.lc, init.years, iter.max, stop.crit, ncores, print.log=TRUE, plot=TRUE, writeRaster=TRUE, korr_iter=1)
 
argument | description 
----- | ----
lc | categorical RasterLayer of the initial Land Use/Cover Classes  
suit | either a RasterStack or a list of RasterStacks(for each year) of the probabilities for the modeled land use/cover classes resulting from the statistical modelling. The data type should be Float (FLT4S). The names of the layers should correspond to the landuse/cover classes as follows: "lc1", "lc2", "lc3",..  
spatial | either a RasterLayer or a list of RasterLayers(for each year) of locations where no land sue change is allowed (i.e. Protected Areas) containing the values 0 for c areas where conversions are allowed and 1 for areas where conversions are not allowed
demand | matrix specifying the amount of pixel for each land use/cover class in the subsequent modelling steps. Columns are land use/cover classes, number of rows equal the number of modelling steps. Values should be integer.
elas | vector containing values referring to the conversion elasticity of the land use/cover classes. 0: easy to convert, 0.5 : medium to convert, 1: difficult to convert.
traj | matrix describing the trajectories of land use/cover. Rows: initial land use/cover, Columns: following land use/cover. Values define the years of transition, e.g. 0: no transition allowed, 1: transition allowed after first iteration, 10: transition allowed after 10 iterations.
nochange.lc | vector with integer numbers of stable classes, e.g. water.
init.years | factor to set the initial number of years the pixels are under the specific land use/cover at the beginning of the modelling.
stop.crit | vector containing 3 values. the first one to the max deviation of allocated land use/cover to the demand, the second one to the maximum deviation of pixels for the smallest demand class, and the third to the maximum deviation of each demand class.
iter.max | integer number specifying the maximum number of iteration until the allocation of land use/cover is stopped
ncores | integer number specifying the number of cores to us during processing
print.log | TRUE/FALSE if tail of log file is printed during processing
writeRaster | TRUE/FALSE if scenario output raster should be written to the working directory during iteration

## Output: 
list of:    
[[1]] "RasterStack" containing the categorical scenarios of land use/cover allocation for the requested years   
[[2]] "matrix" of all log information (pixel difference, correction value, iterValue)   
[[3]] "matrix" of the final log information of the last allocation for each epoche (pixel difference, correction value, iterValue)   

##Procedure
1. Non-Spatial Domain:  
  the amount of land use/cover change for the land use/cover scenarios (row numbers represent modelling steps) has to be provided externally in form of a matrix. Vaules correspond to Pixel of land use/cover in the respective year. Options to derive possible change of land use/cover are manifold. One of the most simple assumptions of future land use/cover amount could follow a simple extrapolation of trends from known from "historic" datasets or similar informations.

2. Spatial Domain:  
2.1. alucR takes as input data a current land-cover map and a raster stack of land use/cover "suitability" for each land use/cover class for which we want to generate   the scenario. The suitability raster usually result from a statistical analysis, e.g logistic regression (to estimate the probability of a pixel to be under a certain land use/cover), or similar approaches, representing continuous values between 0 an 1.  
  A preliminary land use/cover  scenario is generated assigning the land use/cover class which has the highest class probability from the suitability stack for each pixel to the scenario map. The amount of land use/cover is then compared with the requested demand and the suitability layer are iteratevily weighted until the scenario land use/cover classes ~equals the demand requirements. The stop criterium is currently defined as either deviating not more than 0.3percent from the demanded amount of land use/cover plus, having not more than 1 pixel difference for the smallest land use/cover class, or deviating less than 10 pixels, plus  having not more than 1 pixel difference for the smallest land use/cover class.      
2.2. Elasticities refer to the conversion impediments of one land use/cover to another. Example: It might be easy to convert pasture areas to cropland, but very cost intensive to covert urban areas to pasture.  
2.3. Spatial Restrictions refer to Protected areas or other areas where no change of land use/cover is allowed during the simulation. These areas are masked before the new allocation of land use/cover is calculated and later added to scenario rasters.  
2.4. Trajectories define after how many years a land use/cover class can change to another land use/cover class and, if set to 0 that no-change to that land use/cover class is allowed. Example: to change from secondary forest to primary fores we may define a minimum of 10 years to pass before we consider the ixel primary forest again. Or we may not allow urban to change back to another land use/cover class.  
2.5. No change or stable classes are masked out at the beginning of the function and added at the end. You do not need to provide a demand for these classes. Examples might be "water" or "clouds" (resulting from a land cover classification).  
2.6. Initial years defines the amount of years the current land use is there already (relevant for trajectories 2.4.).
  
## Outlook
The current code provides a snapshot of the development and will further be updated and documented. 

### Reference:  
Verburg PH, Soepboer W, Veldkamp A, Limpiada R, Espaldon V, Mastura, Sharifah S. A. (2002) Modeling the Spatial Dynamics of Regional Land Use: The CLUE-S Model. Environmental Management, vols 30(3):391-405