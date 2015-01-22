library(rgdal)

library(raster)


setwd("C:/Users/geo_flgo/Documents/GitHub/alucR/example_data/small")
lc<- raster("tc08_aggregated.tif" )
lc10 <- raster("tc10_aggregated.tif" )

lcPix08 <- table(getValues(lc))
lcPix10 <- table(getValues(lc10))


diff.year <- (lcPix10 - lcPix08)/2
diff.year

demand_agg <- rbind(lcPix08 + diff.year)
demand_agg
for (i in 2:22){
  demand_agg <- rbind (demand_agg, demand_agg[(i-1),] + diff.year)
}
demand_agg
demand <- round(demand_agg)
demand 
write.csv (demand, "demand_aggregated.csv")
