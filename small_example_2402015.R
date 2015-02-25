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
suit<- stack("lu_suitability_stack.tif"); names(suit)<- c("lc7","lc4","lc3")
spatial<- raster("PAall_aggregated.tif")

#load csv
demand<- read.csv("lu_demand_aggregated.csv", row.names=1)
names(demand)<- c("lc3", "lc4", "lc7")

nochange.lc <- c("lc5","lc6","lc8") # müssen ausmaskiert werden und später wieder eingesetzt
natural.lc <- c("lc1","lc2")

elas<- read.csv("elas.csv",row.names=1)
traj <- read.csv("trajectories.csv",row.names=1)
init.years=5
ncores=detectCores()/2

# sceanrios
scenario <- aluc(lc=lc, suit=suit, natural=c("lc1","lc2"), nochange.lc=c("lc5","lc6","lc8"), spatial=spatial, demand=demand, elas=as.matrix(elas), traj=as.matrix(traj),
                 init.years= 5,  stop.crit=c(0.0003 , 1, 10), iter.max=100, ncores=detectCores()/2, print.log=TRUE, print.plot=TRUE, writeRaster=FALSE)

#view output
scenario.raster <- scenario[[1]] 
plot(scenario.raster,y=1:16, col=c("darkgreen", "lightgreen", "orange", "yellow","white","darkblue" , "red", "grey"))
logfile1 <- scenario[[2]]
View(logfile1)

