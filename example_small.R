#author: florian.gollnow@geo.hu-berlin.de


#example script:
setwd(".../small")

lc<- raster("tc08_aggregated.tif" )
suit<- stack("suitability_stack.tif"); names(suit)<- c("lc1","lc2","lc3","lc4","lc7","lc8")
spatial<- raster("PAall_aggregated.tif")

demand<- read.csv("demand_aggregated.csv", row.names=1)
elas<- read.csv("elas.csv",row.names=1)
traj <- read.csv("trajectories.csv",row.names=1)


scenarios <- aluc(lc=lc, suit=suit, spatial=spatial,demand=as.matrix(demand[c(1:3),]), elas=as.matrix(elas), traj=as.matrix(traj), nochange.lc=c(5,6), init.years=5, ncores=detectCores()/4,iter.max=100, writeRaster=FALSE, korr_iter=1)

scenario.raster <- scenarios[[1]] 

plot(scenario.raster)