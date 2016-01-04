
library(raster)
setwd("C:/Users/geo_flgo/Documents/GitHub/alucR/example_data/big")

dir()


lc08<-raster("tc08.tif")
lc10<-raster("tc10.tif")

lcPix08 <- table(getValues(lc08))
lcPix10 <- table(getValues(lc10))


diff.year <- (lcPix10 - lcPix08)/2
diff.year

demand_lin <- rbind(lcPix08 + diff.year)
demand_lin
for (i in 2:22){
  demand_lin <- rbind (demand_lin, demand_lin[(i-1),] + diff.year)
}
demand_lin
demand <- round(demand_lin)
demand 
write.csv (demand, "demand.csv")
