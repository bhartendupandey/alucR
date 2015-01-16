
library(rgdal)
library(sp)
library(raster)
library(parallel)

setwd("q:/florian/data/studyclue/r_model_2015")
setwd("q:/carbiocial/florian/data/studyclue/r_model_2015")

#tables
demand<- read.csv("demand.csv", header=F)
elas <- read.csv("elasticities.csv", header=F)
trajectories <-  read.csv("trajectories.csv", header=F)
nochange.lc=c(5,6)

init.years <- 5

#raster
lc <- raster("tc08.tif"); lc[Which (lc==0 | lc==9 ) ] <- NA 

suit <- stack("pforestES.tif", "psecvegES.tif", "ppastureES.tif","pcropES.tif", 
              "urban_buffer.tif", "potherES.tif") 
names(suit)<- c("lc1", "lc2", "lc3", "lc4", "lc7", "lc8")

spatial <- raster("pa_all.tif")

#
ncores = detectCores()/4 
iter.max=110
max.epoche =nrow(demand)

#####
#packages for processing
library(rgdal)
library(sp)
library(raster)
library(parallel)
####

scenarios<-aluc(lc=lc, suit=suit, spatial=spatial,demand=demand, elas=elas, traj=trajectories, nochange.lc=nochange.lc, init.years=5, ncores=detectCores()/4,iter.max=100, writeRaster=TRUE)


aluc <- function (lc, suit, spatial, demand, elas, traj, nochange.lc, init.years, ncores, iter.max, print.log=TRUE, plot=TRUE, writeRaster=TRUE){

epoche=1

while (epoche <= nrow(demand)){
#####
#raster to vector  
data_vector <- if (epoche==1) {getValues(lc)
                  }else {
                    new.data_vector
                  }

p_vector <-   if(class(suit)=="RasterStack"){ 
                  getValues(suit)
              }else if (class(suit)=="character"){
                  getValues(stack(get(suit)[epoche]))
                  }             
  
sp.rest_vector <- if(class(spatial)=="RasterLayer"){ 
                      getValues(spatial)
                  }else if (class(spatial)=="character"){
                       getValues(stack(get(spatial)[epoche]))
                        }        

#####
print ("raster to vector conversion done")					
#####
lc_cat <- as.numeric(gsub("lc","",colnames(p_vector)))   # tolower(colnames(p_vector))
layer <- ncol(p_vector)


#mask NA's from p_vector
p_vector [is.na(data_vector), ] <- NA

#####
# only for the first epoche, forthe following epoches defined below:
# years under land use only for the first epoche, forthe following epoches defined below
# previous land use
if (epoche==1){
    trans.years_vector <- rep(init.years, length (data_vector))
    tprop.previous_vector <- data_vector
  }                   
#
#####
print ("start trajectories")
#####
#Trajectories
print("start trajectories")
tprop.previous_landuse=which(!is.na(tprop.previous_vector)) # only id's with lc (no NAs)
tprop.previous_layers=as.vector(na.exclude(unique(tprop.previous_vector)))# unique Land cover
for (i in tprop.previous_landuse) { # loop through all pixels with lc
 lu=tprop.previous_vector[i];
  for (t in tprop.previous_layers) { # loop through all layers of suitability
    if (trans.years_vector[i] < traj[lu,t]) {
      p_vector[i,t]=NA;
    }
  }
}
print("trajectories/matrix - done")  

#####
print ("spatial restrictions")
#####
#Spatial restrictions
sp.rest_index <- which(!is.na(sp.rest_vector));
p_vector[sp.rest_index,] <- NA;

#adjust demand to spatial restrictions
lc.sp.rest <- tabulate(data_vector[sp.rest_index], nbins=max(lc_cat))

#noChange classes
demand.new <- demand[epoche,]- lc.sp.rest
demand.new [nochange.lc] <- 0 # set demand to zero for no change classes

#set no.change classes to NA
nochange_index <- is.element(data_vector, nochange.lc) 
p_vector[nochange_index, ]<- NA

#####
print ("prepare suitability array")
#####
#normalize p_vector
for (i in c(1:layer)) {
  stretch <- 1/max(p_vector[,i],na.rm=TRUE);
  p_vector[,i] <- p_vector[,i]*stretch;
}

# add ELAS auf P for every layer fehler in loop

for (i in 1:length(lc_cat)) {
  elas_index <- which(data_vector==lc_cat[i])
  temp <- p_vector[elas_index, i] + as.numeric(elas[lc_cat[i]])
  p_vector[elas_index,i] <- temp
}

#####
# set initial values for ITER 
if (epoche==1) {
  iter <- rep(0,max(lc_cat));
} else {
  iter <- iter.hist[nrow(iter.hist),]; #from the last epoche  
}

#####
# prepare log data frames
if (epoche==1){
  log1 <- c()
  log2 <- c()
  names.log <- c("u", paste("pix.diff_lc",lc_cat, sep=""), paste("korr_lc",lc_cat, sep=""), paste("iter_lc",lc_cat, sep=""))
}

iter.hist <- c(iter)
diff.pix.hist <- c()
diff.p.hist <- c()
change.p.hist <- c()
korr.hist <- c()

#####
#start allocation iteration
#####
print (paste("start iteration for epoche:", epoche, sep=""))
#####
u=1
#multicore definition
cl <- makeCluster(getOption("cl.cores", ncores));

repeat { 
  # calc P2 values with ITER 
  if (u==1 & epoche==1){
    p2_vector<- p_vector
  }else{
    for (i in c(1:layer)) {
      p2_vector[,i] <- p_vector[,i]+as.numeric(iter[lc_cat[i]]);
    }
  }
  
  # find winner
  # if all layer have NA, dann NA, sonst which.max(P2)
  tprop_vector <- lc_cat[parRapply(cl,p2_vector,FUN=function(w) ifelse(all(is.na(w)),NA,which.max(w)))];
  #tprop_temp <- parRapply(cl,p2_vector,FUN=function(w) ifelse(all(is.na(w)),NA,which.max(w)))
  #tprop_vector<- ifelse(is.na(tprop_temp), NA, lc_cat[tprop_temp])
  n <- tabulate(tprop_vector,nbins=max(lc_cat))
  
  diff.pix <- n-demand.new
  diff.perc <- diff.pix/demand.new
  #save history
  diff.pix.hist <- rbind(diff.pix.hist, diff.pix)
  diff.p.hist <- rbind(diff.p.hist, diff.perc)
  
  
  #korr
  if(u==1){ # initializing Korr 
    korr= -1*sign(diff.perc)*1/100
    
  }else{ # modifying Korr based on prior results
    change.perc <- abs((diff.pix.hist[u-1,]-diff.pix.hist[u,])/(diff.pix.hist[u-1,])*100) 
    proportion <- abs(diff.pix.hist[u,])/ colSums(abs(diff.pix.hist[c(u,u-1),]),na.rm=TRUE)
    better<- abs(diff.pix.hist[u,])< abs(diff.pix.hist[u-1,])
    
    #korr <- as.vector(ifelse(diff.pix.hist[u,]== 0 , 0,
     #                        ifelse(diff.pix.hist[u,]!=0 & korr.hist[u-1,]==0,-1*sign(diff.perc)*1/100,
      #                              ifelse(better  ==FALSE & sign(diff.pix.hist[u,])  ==sign(diff.pix.hist[u-1,]), korr.hist[u-1,]*2,
       #                                    ifelse(abs(diff.perc)<= 0.001 & abs(diff.pix.hist[u,])<= 20 & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), (korr.hist[u-1,]/2)+(korr.hist[u-1,]/2)*proportion,
        #                                          ifelse(change.perc< 20& sign(diff.perc)==sign(diff.p.hist[u-1,]),korr.hist[u-1,]*2,
         #                                                ifelse(change.perc< 40& change.perc >= 20 & sign(diff.perc)==sign(diff.p.hist[u-1,]), korr.hist[u-1,]+(korr.hist[u-1,]*proportion),
          #                                                      ifelse(sign(diff.pix.hist[u,])	!=sign(diff.pix.hist[u-1,]), -1* korr.hist[u-1,]*proportion,
           #                                                            korr.hist[u-1,]))))))), mode="numeric") 
    
   
    korr <- as.vector(ifelse(diff.pix.hist[u,]== 0 , 0,
                             ifelse(diff.pix.hist[u,]!=0 & korr.hist[u-1,]==0,-1*sign(diff.perc)*1/100,
                                    ifelse(better==FALSE & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), korr.hist[u-1,]*2,
                                           ifelse(abs(diff.perc)<= 0.001 & abs(diff.pix.hist[u,])<= 20 & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), (korr.hist[u-1,]/2)+(korr.hist[u-1,]/2)*proportion,
                                                  ifelse(change.perc< 40& sign(diff.perc)==sign(diff.p.hist[u-1,]), korr.hist[u-1,]+(korr.hist[u-1,]*proportion),
                                                         ifelse(sign(diff.pix.hist[u,])  !=sign(diff.pix.hist[u-1,]), -1* korr.hist[u-1,]*proportion,
                                                                korr.hist[u-1,])))))), mode="numeric") 
    
    
    
    change.p.hist <- rbind(change.p.hist, change.perc)
  }
  
  iter <- iter + korr;
  
  #save history
  korr.hist <- rbind(korr.hist, korr)
  iter.hist <- rbind(iter.hist, iter)

#####  
  
  if(plot==TRUE){
    plot(0,0,xlim = c(2,iter.max),ylim = c(-1.0,1.0),ylab="iter", xlab="iteration", type = "n")
    legend("topright", legend=paste("LC",lc_cat, sep=""), col=rainbow(length(lc_cat)), pch=15)
    for (i in 1:length(lc_cat)){
      lines(c(1:nrow(iter.hist)),iter.hist[,lc_cat[i]],col=rainbow(length(lc_cat))[i],type = 'l', lwd=1.5);
    }
  }
  
  log.tmp <- as.vector(c(u, diff.pix[lc_cat], korr[lc_cat], iter[lc_cat]), mode="numeric")
  names(log.tmp) <- names.log

  log1 <- rbind(log1,log.tmp)
  
  if(print.log==TRUE){
    print(tail(log1))
  }
#####
  #initialize next u sequence
  u=u+1;
  #stop argument iteration 
  if (max(abs(diff.perc),na.rm=TRUE) < 0.001 & abs(diff.pix [7]) <= 1 & max(abs(diff.pix),na.rm=TRUE)<=10)  {
    break;
  }
  if (abs(diff.pix [7]) <= 1 & max(abs(diff.pix),na.rm=TRUE)<=10){
    break;
  }
  # stop argument iteration if ITERmax reached
  if (u > iter.max) { 
    break;
  }
}# next iteration over "u"
#stop cluster
stopCluster(cl);
print("allocation done")
log2 <- rbind(log2,log1[nrow(log1),])
#####
#insert spatial restrictions
tprop_vector[sp.rest_index] <- data_vector[sp.rest_index];

#stable classes
#nochange_index <- match(nochange.lc, data_vector)
tprop_vector[nochange_index] <- data_vector[nochange_index]


## compare this allocation for transistion years, inc if changed, reset to 1 if change
transition_years_vector <- ifelse(tprop_vector==tprop.previous_vector, trans.years_vector + 1, 1);

# save as raster object 
scenario <- lc
scenario <- setValues(scenario, tprop_vector)

if (writeRaster==TRUE){writeRaster(scenario, paste("scenario", epoche, ".tif", sep=""))}

assign(paste("scenario", epoche, sep=""), scenario)

  
# now set previous to this epoche and start next
tprop_previous_vector <- tprop_vector;

print("epoche done")
} # end of epoche loop 

return(list(stack (mget (paste("scenario", rep(1:nrow(demand)),sep=""))), log1, log2))
}
  
  
  