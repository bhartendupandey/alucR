# file:     example_allocate_v2.R
#
# coder:    
# florian.gollnow@geo.hu-berlin.de
# 

library(rgdal)
library(sp)
library(raster)
library(parallel)

#setwd("q:/florian/data/studyclue/r_model_2015")
setwd("q:/carbiocial/florian/data/studyclue/r_model_2015")

#raster
lc <- raster("tc08.tif"); lc[Which (lc==0 | lc==9 ) ] <- NA 
suit <- stack("pforestES.tif", "psecvegES.tif", "ppastureES.tif","pcropES.tif", 
              "urban_buffer.tif", "potherES.tif") 
names(suit)<- c("lc1", "lc2", "lc3", "lc4", "lc7", "lc8") #naming scheme
spatial <- raster("pa_all.tif")

#tables & vectors
demand<- read.csv("demand.csv", header=F)[c(1:5),]
elas <- read.csv("elasticities.csv", header=F)
trajectories <-  read.csv("trajectories.csv", header=F)
nochange.lc=c(5,6)
#



# run land use model:
scenarios<-aluc(lc=lc, suit=suit, spatial=spatial,demand=demand, elas=elas, traj=trajectories, nochange.lc=nochange.lc, init.years=5, ncores=detectCores()/4,iter.max=100, writeRaster=TRUE)


scenario_raster <- scenarios[[1]]


#plot the scenario_raster 
plot(scenario_raster)


