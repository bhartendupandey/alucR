# file:     allocate_v2.R
#
# coder:    
# florian.gollnow@geo.hu-berlin.de
# moenkemt@geo.hu-berlin.de



#alucR - allocation of land use change 
#---
## alucR
#alucR - Project is a first step to implement a Land Use Change Model in R (http://www.r-project.org). 
#We have been following the basic framework provided by Verburg et al. (2002). The code uses basic R-language and packages and is fully documented. This makes
#it possible to easily adapt the code to the users specific needs. 


## Function definition:
# aluc(lc, suit, spatial, demand, elas, traj, nochange.lc, init.years, ncores, iter.max, print.log=TRUE, plot=TRUE, writeRaster=TRUE)

 
#argument | description 
#----- | ----- 
#lc | categorical RasterLayer of the initial Land Use Classes  
#suit | either a RasterStack or a list of RasterStacks(for each year) of the probabilities for the modeled land use classes resulting from the statistical modelling. The datatype should be Float (FLT4S). The names of the layers should correspond to the landuse/cover classes as follows: "lc1", "lc2", "lc3",..  
#spatial | either a RasterLayer or a list of RasterLayers(for each year) of locations where no land sue change is allowed (i.e. Protected Areas) containing the values 0 for c areas where conversions are allowed and 1 for areas where conversions are not allowed
#demand | matrix specifying the amount of pixel for each land use class in the subsequent modelling steps. Columns are land use classes, number of rows equal the number of modelling steps. Values should be integer.
#elas | vector containing values referring to the conversion elasticity of the land use/cover classes. 0: easy to convert, 0.5 : medium to convert, 1: difficult to convert.
#traj | matrix describing the trajectories of land use. Rows: initial land use/cover, Columns: following land use/cover. Values define the years of transition, e.g. 0: no transition allowed, 1: transition allowed after first iteration, 10: transition allowed after 10 iterations.
#nochange.lc | vector with integer numbers of stable classes, e.g. water.
#init.years | factor to set the initial number of years the pixels are under the specific land use at the beginning of the modelling.
#iter.max | integer number specifying the maximum number of iteration until the allocation of land use is stopped
#ncores | integer number specifying the number of cores to us during processing
#print.log | TRUE/FALSE if tail of log file is printed during processing
#writeRaster | TRUE/FALSE if scenario output raster should be written to the working directory during iteration

## Output: 
#list of    
#[[1]] RasterStack containing the categorical scenarios of land use allocation for the requested years   
#[[2]] matrix of all log information   
#[[3]] matrix of final log information of the final allocation for each epoche.   


#####
#packages dependencies
library(rgdal)
library(sp)
library(raster)
library(parallel)
####

aluc <- function (lc, suit, spatial, demand, elas=rep(0, max(lc_cat)), traj=matrix(data=1, ncol=max(lc_cat), nrow=max(lc_cat)), nochange.lc=c(), init.years= 5,  stop.crit=c(0.003 , 1, 10),iter.max=100, ncores=detectCores(), print.log=TRUE, plot=TRUE, writeRaster=TRUE, korr_iter=1){

epoche=1

while (epoche <= nrow(demand)){
#####
#raster to vector  
data_vector <- if (epoche==1) {
                  getValues(lc)
                }else {
                    getValues (new.data)
                  }

p_vector <-   if(class(suit)=="RasterStack" | class(suit)=="RasterBrick"){ 
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
if (all(complete.cases(data_vector)==FALSE)){
p_vector [is.na(data_vector), ] <- NA
}

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

#ptm <- proc.time()
for (i in 1:layer){ 
  traj_ind <- which (traj[lc_cat[i],lc_cat] != 1) # identify classes with restricted trajectories
  cat_index <- which(tprop.previous_vector==lc_cat[i]) # index classes
  for (a in 1:length(traj_ind)){
    p_vector[cat_index, traj_ind[a]]<- ifelse (trans.years_vector[cat_index] < traj[traj_ind[a], lc_cat[i]], NA, p_vector[cat_index, traj_ind[a]])
    cat(i)
  }
}
#ptm2<- proc.time() - ptm
#print(ptm2)

print("trajectories/matrix - done")  

#####
print ("spatial restrictions")
#####
#Spatial restrictions
sp.rest_index <- which(!is.na(sp.rest_vector));
p_vector[sp.rest_index,] <- NA;

#adjust demand to spatial restrictions
lc.sp.rest <- tabulate(data_vector[sp.rest_index], nbins=max(lc_cat))

#noChange classes - adjust demand 
demand.new <- demand[epoche,]- lc.sp.rest
demand.new [nochange.lc] <- NA # set demand to zero for no change classes
    
min.demand <- which.min(demand.new)
 
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
#if (epoche==1) {
  iter <- rep(0,max(lc_cat));
#} else {
#  iter <- iter.hist[nrow(iter.hist),]; #from the last epoche  
#}

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
 # if (u==1 & epoche==1){
 if (u==1){
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
 
 if (korr_iter ==1){
    korr <- as.vector(ifelse(diff.pix.hist[u,]== 0 , 0,
                             ifelse(diff.pix.hist[u,]!=0 & korr.hist[u-1,]==0,-1*sign(diff.perc)*1/sample(50:150, 1),
                                    ifelse(better==FALSE & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), korr.hist[u-1,]*2,
                                           ifelse(abs(diff.perc)<= 0.001 & abs(diff.pix.hist[u,])<= 20 & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), (korr.hist[u-1,]/2)+(korr.hist[u-1,]/2)*proportion,
                                                  ifelse(change.perc< 20& sign(diff.perc)==sign(diff.p.hist[u-1,]),korr.hist[u-1,]*2,
                                                         ifelse(change.perc< 40& change.perc >= 20 & sign(diff.perc)==sign(diff.p.hist[u-1,]), korr.hist[u-1,]+(korr.hist[u-1,]*proportion),
                                                                ifelse(sign(diff.pix.hist[u,])  !=sign(diff.pix.hist[u-1,]), -1* korr.hist[u-1,]*proportion,
                                                                       korr.hist[u-1,]))))))), mode="numeric") 
    }else if (korr_iter == 2){
	korr <- as.vector(ifelse(diff.pix.hist[u,]== 0 , 0,
                             ifelse(diff.pix.hist[u,]!=0 & korr.hist[u-1,]==0,-1*sign(diff.perc)*1/sample(50:150, 1),
                                    ifelse(better==FALSE & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), korr.hist[u-1,]*2,
                                           ifelse(abs(diff.perc)<= 0.001 & abs(diff.pix.hist[u,])<= 20 & sign(diff.pix.hist[u,])==sign(diff.pix.hist[u-1,]), (korr.hist[u-1,]/2)+(korr.hist[u-1,]/2)*proportion,
                                                       ifelse(change.perc < 40 & sign(diff.perc)==sign(diff.p.hist[u-1,]), korr.hist[u-1,]+(korr.hist[u-1,]*proportion),
                                                                ifelse(sign(diff.pix.hist[u,])  !=sign(diff.pix.hist[u-1,]), -1* korr.hist[u-1,]*proportion,
                                                                       korr.hist[u-1,])))))), mode="numeric") 
	}
     
    change.p.hist <- rbind(change.p.hist, change.perc)
  }
  
  iter <- iter + korr ;
  iter <- as.numeric (ifelse(iter <=-1, -1, ifelse(iter>=1,1, iter)))
  
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
    #print(tail(log1))
    print(log.tmp)
  }
#####
  #initialize next u sequence
  u=u+1;
  #stop argument iteration 
  
  if (max(abs(diff.perc),na.rm=TRUE) < stop.crit[1] & abs(diff.pix [min.demand]) <= stop.crit[2] )  {
    break;
  }
  if (abs(diff.pix [min.demand]) <= stop.crit[2])& max(abs(diff.pix),na.rm=TRUE)<=stop.crit[3])){
    break;
  }
  # stop argument iteration if ITERmax reached and take the ITER with the minimum deviation from the demand from all iterations 
    if (u > iter.max) {
        current.log <- log1[(nrow(log1)-iter.max+1):nrow(log1),]
        iterfinal_index <-  which.min(rowSums(abs(current.log[,2:(layer+1)])))
        iterfinal <- current.log[iterfinal_index, (2*layer+2):(3*layer+1)]
        print(iterfinal)
        for (i in 1:layer){ 
			p2_vector[,i] <- p_vector[,i]+ as.numeric(iterfinal[i]);
        }	    
        log.tmp <- as.vector(c(I(iter.max+1), current.log [iterfinal_index,2:(3*layer+1)]), mode="numeric")
        log1 <- rbind(log1, log.tmp)
        print("break")
        print(log.tmp )
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
new.data <- lc
new.data <- setValues(new.data, tprop_vector)

if (writeRaster==TRUE){writeRaster(new.data, paste("scenario", epoche, ".tif", sep=""), overwrite=TRUE)}

assign(paste("scenario", epoche, sep=""), new.data)

  
# now set previous to this epoche and start next
tprop_previous_vector <- tprop_vector;

print("epoche done")
#initialize new epoche
epoche <- epoche+1
} # end of epoche loop 

return(list(stack (mget (paste("scenario", rep(1:nrow(demand)),sep=""))), log1, log2))
}
