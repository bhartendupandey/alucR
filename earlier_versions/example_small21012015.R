#author: florian.gollnow@geo.hu-berlin.de

# Landuse/cover classes 
# LC1 Forest
# LC2 Secondary Vegetation
# LC3 Pasture
# LC4 Cropland
# LC5 Clouds
# LC6 Water
# LC7 Urban
# LC8 Other
# The land use/cover data is a reclassified TerraClass product of INPE http://www.inpe.br/cra/projetos_pesquisas/terraclass2010.php    

#example script:
#set workingdirectory
setwd("C:/Users/geo_flgo/Documents/GitHub/alucR/example_data/small")

#load packages
library(rgdal)
library(sp)
library(raster)
library(parallel)

#load raster
lc<- raster("tc08_aggregated.tif" )
suit<- stack("suitability_stack.tif"); names(suit)<- c("lc1","lc2","lc3","lc4","lc7","lc8")
spatial<- raster("PAall_aggregated.tif")

#load csv
demand<- read.csv("demand_aggregated.csv", row.names=1)
elas<- read.csv("elas.csv",row.names=1)
traj <- read.csv("trajectories.csv",row.names=1)

# sceanrios
scenarios <- aluc(lc=lc, suit=suit, spatial=spatial,demand=as.matrix(demand), elas=as.matrix(elas), traj=as.matrix(traj), nochange.lc=c(5,6), init.years=5, ncores=detectCores()/2,iter.max=100, writeRaster=FALSE, korr_iter=1)

#view output
scenario.raster <- scenarios[[1]] 
plot(scenario.raster,y=1:16, col=c("darkgreen", "lightgreen", "orange", "yellow","white","darkblue" , "red", "grey"))
log1 <- scenarios[[2]]
View(log1)
log2 <- scenarios[[3]]
View (log2)
