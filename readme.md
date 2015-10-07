alucR - allocation of land use change 
---

alucR - Project is a first step to implement a Land Use Change Model in R (http://www.r-project.org). We have been following the basic framework provided by Verburg et al. (2002). Land use is hereby spatially allocated based on its location based competitive advatages following the suitability of a certain cell for the specific land use. Suitability might be assessd using statistical methods (for example) or other modelling techniques. The amount of land use to be allocated has to be estimated for the total study area and provided as numbers of pixels. Natural land cover and possible succession stages are modelled based on the temporal trajectories of succession stages defined before in the trajectories matrix. 
The code uses basic R-language and packages and is fully documented. This makes it possible to easily adapt the code to the users specific needs. 

 TO USE THE FUNCTION run the code alucR_function.R in your R-console.   
 Dependencies: raster package; parallel package; rgdal package; sp package

# Function definition:
 aluc (lc, suit, natural, nochange.lc,spatial, demand, elas,traj,init.years, method, stop.crit,iter.max,ncores, print.log,print.plot,writeRaster)

argument | description 
----- | ----- 
lc | categorical RasterLayer of the initial Land Use/Cover Classes  
suit | either a RasterStack or a list of RasterStacks(for each year) of the suitabilities for land cover classes (ordered by preferences) resulting from he statistical modelling. The data type should be Float (FLT4S). The names of the layers should correspond to the landuse classes as follows: "lc7", "lc4", "lc3",..  
natural | character string defining land cover classes refering to natural vegetation ordered by succession states. example: c("lc1", "lc2")
nochange.lc |  character string defining land cover/use classes without suitability layer which are expected to stay unchanged (for example: water). 
spatial | either a RasterLayer or a list of RasterLayers(for each year) of locations where no land use change is allowed (i.e. Protected Areas) containing NA for areas where conversions are allowed and 1 for areas where conversions are not allowed
demand | matrix specifying the amount of pixel for each land use class in the subsequent modelling steps. Columns refer to the land use classes for which there is a suitability layer, number of rows equal the number of modelling steps. Values should be integer.
elas | vector containing values referring to the conversion elasticity of the land use/cover classes. There must be specified for all classes in the land cover product. 0: easy to convert, 0.5 : medium to convert, 1: difficult to convert.
traj | matrix describing the trajectories of land use/cover. Rows: initial land use/cover (1 to n), Columns: following land use/cover (1 to n). Values define the years of transition, e.g. 0: no transition allowed, 1: transition allowed after first iteration, 10: transition allowed after 10 iterations. must be specified for all land_cover classes.
init.years | numeric value or RasterLayer to set the initial number of years the pixels are under the specific land use/cover at the beginning of the modelling.
method | either "competitive" or "hierarchical"
stop.crit | only applies for method="competitive": vector containing 3 values. the first one defines the maximum deviation of allocated land use/cover to the demand in percent, the second one the maximum deviation of pixels for the smallest demand class, and the third defines the maximum deviation of each demand class in pixel.
iter.max | only applies for method="competitive":integer number specifying the maximum number of iteration until the allocation of land use/cover is stopped (in that case the best out of the available allocation is returned)
ncores | only applies for method="competitive":integer number specifying the number of cores to use during processing
print.log | TRUE/FALSE if tail of log file is printed during processing
print.plot | TRUE/FALSE if iter and the final raster are plotted during model execution
writeRaster | TRUE/FALSE if scenario output raster should be written to the working directory during iteration

## Output: 
list of:    
[[1]] RasterStack containing the categorical scenarios of land use allocation for the requested years   
[[2]] matrix of all log information    

##Procedure
1. Non-Spatial Domain - amount of land use to allocate:  
  the amount of land use to be allocated for the land use/cover scenarios (row numbers represent modelling steps) has to be provided externally in form of a matrix. Vaules correspond to Pixel of land use in the respective sceanrio iteration. Options to derive possible scenarios of land use amounts are manifold. One of the most simple assumptions of future land use could follow a simple extrapolation of trends known from "historic" datasets or similar informations.

2. Spatial Domain:  
2.1. alucR takes as input data a current land-cover map and a raster stack of land use "suitability" layer for each land use class we want to generate scenarios for. The suitability raster usually result from a statistical analysis, e.g logistic regression (to estimate the probability of a pixel to be under a certain land use/cover), or similar approaches, representing continuous values between 0 an 1 (Make sure not to have to many similar suitabilities for too many pixels).      
Natural land cover is modeled as a function of land use demands (Natural vegetation = Total Area - Land use - Stable Classe (i.e. water)) and natural success stages defined as temporal trajectores.  
Method: Competitive allocation: A preliminary land use/cover  scenario is generated assigning the land use class which has the highest class probability from the suitability stack for each pixel to the scenario map. The amount of land use/cover is then compared with the requested demand and the suitability layer are iteratively weighted until the scenario land use classes ~equals the demand requirements. The stop criterium is currently defined as either deviating not more than 0.03% from the amount specified in the demand matrix, plus having not more than 1 pixel difference for the smallest land use/cover class, or deviating less than 10 pixels. If this criteria are not meet the iteration stops after the defined number of iteartion (default=100). In this case the map with the scenario map with smalles summed deviation from the requested demand will be saved.       

 Method hirarchical allocation: land use is allocated following the hirachical prefernces. e.g. urban first, followed by crop, followed by pasture. hirarchical allocation is much faster than the competitive allocation.   
2.2. Elasticities refer to the conversion impediments of one land use/cover to another. To keep stable patterns of land use, elasticities should be high, while low elasticities make a class more the land use pixel more likely to change. Relative differences between the land use classes schould also be considered: for example - It might be relatively easy to convert pasture areas to cropland, but very cost intensive to covert urban back to another land use. Elasticities  need to be defined for each landcover class (also stable classes as 0).   
2.3. Spatial Restrictions refer to protected areas or other areas where no change of land use/cover is allowed during the simulation. These areas are masked before the new allocation of land use/cover is calculated and later added to the scenario rasters. Succession stages of natural vegetation is the only possible change within these areas. Spatial restriction don't need to saty the same during the full simulation period, but can vary from simulation step to the next.   
2.4. Trajectories define after how many years a land use/cover class can change to another land use/cover class. If set to 0no-change to that land use/cover class is allowed. Example: to change from secondary forest to primary fores we may define a minimum of 10 years to pass before we consider the ixel primary forest again. Or we may not allow urban to change back to another land use/cover class.        
2.5. No change or stable classes are masked out at the beginning of the function and added at the end. You do not need to provide a demand for these classes. Examples might be "water" or "clouds" (resulting from a land cover classification).      
2.6. Initial years defines the amount of years the current land use is there already (relevant for trajectories 2.4.).     
3. Structure of code:

* 1. Pre-processing   
  1.1 Read data   
  1.2 Pseudo natural vegetation layer   
  1.3 Descriptive variables   
  1.4 Spatial Restrictions   
	  1.4.1 NA values in land cover dataset   
	  1.4.2 Defined stable (no change) classes   
	  1.4.3 Defined spatial restrictions   
  1.5 Adjusting demand of land use   
	  1.5.1 Add natural land cover demand   
	  1.5.2 Adjust for spatial restrictions (4.3)   
  1.6 Trajectories of land use change   
	  1.6.1 Trajectories for land use classes   
	  1.6.2 Trajectories for natural vegetation   
  1.7 Add elasticities   
  1.8 Combine land use suitability and natural vegetation layer   
* 2. Allocation sub module    
	2.1 Definition of allocation function   
	2.2 Function - allocation module   
		2.2.1 Initiate & start iteration   
		2.2.2 Read data & add iter   
		2.2.3 Evaluate competitive advantages   
		2.2.4 Compare amount of received classes to demand    
		2.2.5 Adjust iter for next iteration (starting 2.2.2)   
		2.2.6 Stop iteration when amount of classes ~ demand 	   
		2.2.7 return allocation vector and logfile   
	2.3 Run allocation module (defined above)   
* 3. Post-processing   
	3.1 Read results from allocation module (2.2)   
	3.2 Spatial restrictions    
		3.2.1 Add stable classes (1.4.2)   
		3.2.2 Spatial restrictions (1.4.3)   
	3.3 Natural vegetation and succession   
		3.3.1 Reclassify all natural to pseudo natural classes   
		3.3.2 Reclassify pseudo natural class based on trajectory and succession order   
	3.4 Saving results and preparing next epoche   
		3.4.1 Updating transition years vector
		3.4.2 Save final allocation result as raster   
* 4. Return results and logfile   
################################      
  
## Outlook
The current code provides a running snapshot of the development and will further be updated and documented. 

##NEXT Steps:
1. stick with rasters instead of convering everything to vectors
2. stop land use allocation if demand cannot be allocated 
3. ...

### Reference:  
Verburg PH, Soepboer W, Veldkamp A, Limpiada R, Espaldon V, Mastura, Sharifah S. A. (2002) Modeling the Spatial Dynamics of Regional Land Use: The CLUE-S Model. Environmental Management, vols 30(3):391-405