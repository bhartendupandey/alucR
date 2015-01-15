alucR - allocation of land use change 
---
# alucR
alucR - Project is a first step to implement a Land Use Change Model in R (http://www.r-project.org). 
We have been following the basic framework provided by Verburg et al. (2002). The code uses basic R-language and packages and is fully documented. This makes
it possible to easily adapt the code to the users specific needs. 
The large community of R-users might additional help to implement improvements in performance and flexibility to the current procedures.


## Function definition:
 aluc(lc, suit, spatial, demand, elas, traj, nochange.lc, init.years, ncores, iter.max, print.log=TRUE, plot=TRUE, writeRaster=TRUE)

 
argument | description 
----- | ----- 
lc | categorical RasterLayer of the initial Land Use Classes  
suit | either a RasterStack or a list of RasterStacks(for each year) of the probabilities for the modeled land use classes resulting from the statistical modelling. The datatype should be Float (FLT4S). The names of the layers should correspond to the landuse/cover classes as follows: "lc1", "lc2", "lc3",..  
spatial | either a RasterLayer or a list of RasterLayers(for each year) of locations where no land sue change is allowed (i.e. Protected Areas) containing the values 0 for c areas where conversions are allowed and 1 for areas where conversions are not allowed
demand | matrix specifying the amount of pixel for each land use class in the subsequent modelling steps. Columns are land use classes, number of rows equal the number of modelling steps. Values should be integer.
elas | vector containing values referring to the conversion elasticity of the land use/cover classes. 0: easy to convert, 0.5 : medium to convert, 1: difficult to convert.
traj | matrix describing the trajectories of land use. Rows: initial land use/cover, Columns: following land use/cover. Values define the years of transition, e.g. 0: no transition allowed, 1: transition allowed after first iteration, 10: transition allowed after 10 iterations.
nochange.lc | vector with integer numbers of stable classes, e.g. water.
init.years | factor to set the initial number of years the pixels are under the specific land use at the beginning of the modelling.
iter.max | integer number specifying the maximum number of iteration until the allocation of land use is stopped
ncores | integer number specifying the number of cores to us during processing
print.log | TRUE/FALSE if tail of log file is printed during processing
writeRaster | TRUE/FALSE if scenario output raster should be writen to the working directory during iteration

## Output: 
list of    
[[1]] RasterStack containing the categorical scenarios of land use allocation for the requested years   
[[2]] matrix of all log information   
[[3]] matrix of final log information of the final allocation for each epoche.   



## Outlook
The current code provides a snapshot of the development and will further be updated and documented. 


 

Reference:
Verburg PH, Soepboer W, Veldkamp A, Limpiada R, Espaldon V, Mastura, Sharifah S. A. (2002) Modeling the Spatial Dynamics of Regional Land Use: The CLUE-S Model. Environmental Management, vols 30(3):391-405