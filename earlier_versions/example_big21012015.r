# example_big.R
#
# florian.gollnow@geo.hu-berlin.de
# 
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

# load packaged
library(rgdal)
library(sp)
library(raster)
library(parallel)

# setwd("q:/florian/data/studyclue/r_model_2015")
setwd("q:/carbiocial/florian/data/studyclue/r_model_2015")

# raster
lc <- raster("tc08.tif")
suit <- stack("pforestES.tif", "psecvegES.tif", "ppastureES.tif","pcropES.tif", 
              "urban_buffer.tif", "potherES.tif") 
names(suit)<- c("lc1", "lc2", "lc3", "lc4", "lc7", "lc8") #naming scheme
spatial <- raster("pa_all.tif")

# tables & vectors
demand<- read.csv("demand.csv", header=F)
elas <- read.csv("elasticities.csv", header=F)
trajectories <-  read.csv("trajectories.csv", header=F)
nochange.lc=c(5,6)


# run land use model:
scenarios<-aluc(lc=lc, suit=suit, spatial=spatial,demand=as.matrix(demand), elas=as.matrix(elas), traj=as.matrix(trajectories), nochange.lc=nochange.lc, init.years=5, ncores=round(detectCores()/4),iter.max=100, writeRaster=FALSE)

# Output
scenario_raster <- scenarios[[1]]
plot(scenario_raster, col=c("darkgreen", "lightgreen", "orange", "yellow","white","darkblue" , "red", "grey"))
log1 <- scenarios[[2]]
View(log1)
log2 <- scenarios[[3]]
View (log2)


