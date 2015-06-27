# file:     alucR_function.R
#
# To use this function just source (run) the script in the R-console. 
#
# coder:    
# florian.gollnow@geo.hu-berlin.de
# moenkemt@geo.hu-berlin.de



#alucR - allocation of land use change 
#---
## alucR
# alucR - Project is a first step to implement a Land Use Change Model in R (http://www.r-project.org). 
# We have been following the basic framework provided by Verburg et al. (2002). Land use is modelled separately to natural vegetation. While land use classes 
# have certain suitability at different locations (depending for example on slope, soil, precipitation or accessibility (road network)) natural vegetation stages are modelled 
# as steps of succession defined in the trajectories (traj) file.  
# The code uses basic R-language and packages and is fully documented. This makes
# it possible to easily adapt the code to the users specific needs. 


## Function definition:
# aluc (lc, suit, natural, nochange.lc,spatial, demand, elas,traj,init.years,method ,stop.crit,iter.max,ncores, print.log,print.plot,writeRaster)

#argument | description 
#----- | ----- 
#lc | categorical RasterLayer of the initial Land Use/Cover Classes  
#suit | either a RasterStack or a list of RasterStacks(for each year) of the suitability for land cover classes (ordered by preferences) resulting from the statistical modelling. The data type should be Float (FLT4S). The names of the layers should correspond to the landuse classes as follows: "lc7", "lc4", "lc3",..  
#natural | character string defining land cover classes referring to natural vegetation ordered by succession states. example: c("lc1", "lc2")
#nochange.lc |  character string defining land cover/use classes without suitability layer which are expected to stay unchanged (for example: water). 
#spatial | either a RasterLayer or a list of RasterLayers(for each year) of locations where no land use change is allowed (i.e. Protected Areas) containing NA for areas where conversions are allowed and 1 for areas where conversions are not allowed
#demand | matrix specifying the amount of pixel for each land use class in the subsequent modelling steps. Columns refer to the land use classes for which there is a suitability layer, number of rows equal the number of modelling steps. Values should be integer.
#elas | vector containing values referring to the conversion elasticity of the land use/cover classes. There must be specified for all classes in the land cover product. 0: easy to convert, 0.5 : medium to convert, 1: difficult to convert.
#traj | matrix describing the trajectories of land use/cover. Rows: initial land use/cover (1 to n), Columns: following land use/cover (1 to n). Values define the years of transition, e.g. 0: no transition allowed, 1: transition allowed after first iteration, 10: transition allowed after 10 iterations. must be specified for all land_cover classes.
#init.years | numeric value or RasterLayer to set the initial number of years the pixels are under the specific land use/cover at the beginning of the modelling.
#method | either "competitive" or "hierarchical"
#stop.crit | only applies for method="competitive": vector containing 3 values. the first one defines the maximum deviation of allocated land use/cover to the demand in percent, the second one the maximum deviation of pixels for the smallest demand class, and the third defines the maximum deviation of each demand class in pixel.
#iter.max | only applies for method="competitive":integer number specifying the maximum number of iteration until the allocation of land use/cover is stopped (in that case the best out of the available allocation is returned)
#ncores | only applies for method="competitive":integer number specifying the number of cores to use during processing
#print.log | TRUE/FALSE if tail of log file is printed during processing
#print.plot | TRUE/FALSE if iter and the final raster are plotted during model execution
#writeRaster | TRUE/FALSE if scenario output raster should be written to the working directory during iteration

## Output: 
#list of    
#[[1]] RasterStack containing the categorical scenarios of land use allocation for the requested years   
#[[2]] matrix of all log information   

######
# Function
######
aluc<-function(  lc, 						
                 suit, 						
                 natural,				
                 nochange.lc=c(),			
                 spatial=c(), 
                 demand=c(), 
                 elas=rep(0, max(lc_unique)), 
                 traj=matrix(data=1, ncol=max(lc_unique, nrow=max(lc_unique))), 
                 init.years= 5,  
                 method = "competitive",
                 stop.crit=c(0.0003 , 1, 10),
                 iter.max=100, 
                 ncores=detectCores()/2, 
                 print.log=TRUE, 
                 print.plot=FALSE, 
                 writeRaster=FALSE) {
 
  ###### Model Structure ############
  #1.Pre-processing
  #  1.1 Read data
  #	1.2 Pseudo natural vegetation layer
  #	1.3 Descriptive variables
  #	1.4 Spatial Restrictions
  #		1.4.1 NA values in land cover dataset
  #		1.4.2 Defined stable (no change) classes
  #		1.4.3 Defined spatial restrictions
  #	1.5 Adjusting demand of land use
  #		1.5.1 Add natural land cover demand
  #		1.5.2 Adjust for spatial restrictions(4.3)
  #	1.6 Trajectories of land use change
  #		1.6.1 Trajectories for land use classes
  #		1.6.2 Trajectories for natural vegetation
  #	1.7 Add elasticities
  #	1.8 Combine land use suitability and natural vegetation layer
  #2. Allocation sub module 
  #	2.1 Definition of allocation function
  #	2.2 Function - allocation module (moved to the end of the script)
  #		2.2.1 Initiate & start iteration
  #		2.2.2 Read data & add iter
  #		2.2.3 Evaluate competitive advantages
  #		2.2.4 Compare amount of received classes to demand 
  #		2.2.5 Adjust iter for next iteration (starting 2.2.2)
  #		2.2.6 Stop iteration when amount of classes ~ demand 	
  #		2.2.7 return allocation vector and logfile
  #	2.3 Run allocation module (defined above)
  #3. Post-processing
  #	3.1 Read results from allocation module 2.2
  #	3.2 Spatial restrictions 
  #		3.2.1 Add stable classes (1.4.2)
  #		3.2.2 Spatial restrictions (1.4.3)
  #	3.3 Natural vegetation and succession
  #		3.3.1 Reclassify all natural to pseudo natural classes
  #		3.3.2 Reclassify pseudo natural class based on trajectory and succession order
  #	3.4 Saving results and preparing next epoche
  #		3.4.1 Updating transition years vector
  #		3.4.2 Save final allocation result as raster
  #4. Return results and logfile
  #################################   
  
##################  
# 1 Preprocessing 
#####  
epoche=1
logfile1 <- c()
while (epoche <= nrow(demand)){
    print(paste("EPOCHE:", epoche , sep=" "))
#####
#  1.1 Read Data
####
    data_vector <- if (epoche==1) {
      getValues(lc)
    }else {
      tprop.previous_vector # getValues (new.data) # change to "tprop.previous_vector" or "tprop_vector" no need to read raster values, since they are stored already
    }
    p_vector <-   if(class(suit)=="RasterStack" | class(suit)=="RasterBrick"){ 
      getValues(suit) # if only one stack is specified
    }else if (class(suit)=="character"){
      getValues(get(suit[epoche])) # in case different stacks for each episode are specified - possibly useful if  for example new roads are build
    }             
    sp.rest_vector <- if(class(spatial)=="RasterLayer"){ 
      getValues(spatial) # spatial restrictions
    }else if (class(spatial)=="character"){
      getValues(get(spatial[epoche])) # in case different stacks for each episode are specified - possibly useful if the protected area network will be expanded during the modelling experiment
    } 
    if (epoche==1){
    trans.years_vector <- if (class(init.years)=="RasterLayer"){
        getValues (init.years)
        } else {
          rep(init.years, length (data_vector))
    }
    if (epoche==1){tprop.previous_vector <- data_vector}
    }
#####
#  1.2 Pseudo natural vegetation layer
#####
    p.natural<- rep(0.5, times=length(data_vector)) # natural vegetation vector 
#####
#  1.3 Descriptive variables
#####
    #land use classes to be modelled
    lu_suit <- as.numeric(gsub("lc","",colnames(p_vector)))   # tolower(colnames(p_vector))
    #no change classes
    no.change <- as.numeric(gsub("lc","",nochange.lc))
    #natural land cover classes (as numeric)
    natural <- as.numeric(gsub("lc","",natural.lc))
    # same as length (lu_suit) number of suitability layers to be modeled
    #lu_suit_l <- ncol(p_vector)
    # pixel for all land use/cover classes
  	lc_pix <- tabulate(data_vector, nbins=max(data_vector, na.rm=TRUE)) 
    # total amount of pixels (excl. NAs)
    lc_n <- sum(lc_pix)
    # unique classes land use/cover classes
    lc_unique <- sort(unique(data_vector))
    #?	# +1 pseudo natural layer for iteration algorithm
    pseudo.N <- max(lc_unique) + 1
    lu.N <- c(lu_suit,  pseudo.N) # class numbers  of all classes to be modelled (incl. pseudo natural class)   
######
#  1.4 Spatial Restrictions   
######
#  	1.4.1 NA values in land cover dataset
#####
    if (all(complete.cases(data_vector))==FALSE){ # skip in case no NAs exits 
      p_vector[is.na(data_vector), ] <- NA
      p.natural[is.na(data_vector)]  <- NA
    } 
#####
#  	1.4.2 Defined stable (no change) classes
#####
    if (length(no.change) > 0 ){
      nochange_index <- is.element(data_vector, no.change) 
      p_vector[nochange_index, ] <- NA
      p.natural[nochange_index]  <- NA
    }  
#####
#  	1.4.3 Defined spatial restrictions
#####
    if (length(sp.rest_vector) > 0 ){ # make sure to 
      sp.rest_index <- which(!is.na(sp.rest_vector));
      p_vector[sp.rest_index,] <- NA;
      p.natural [sp.rest_index] <- NA
    } else { sp.rest_index <- c()}
#####
#  1.5 Adjusting demand of land use
#####
#  	1.5.1 Add natural land cover demand
#####
    nochange.n <- if (length(no.change) > 0 ){
      sum (lc_pix[no.change ])
    }else {c(0)}
    natural.d <- lc_n - nochange.n - sum(demand[epoche,]) 
    print (paste("Area for natural land cover:", natural.d))
#####
#  	1.5.2 Adjust for spatial restrictions(4.3)
#####
    if (length(sp.rest_vector) > 0 ){
    lc.sp.rest <- tabulate(data_vector[sp.rest_index], nbins=max(lc_unique))
    demand.adj <- demand[epoche,] - lc.sp.rest [sort(lu_suit)]
    natural.adj <- natural.d - sum(lc.sp.rest[natural])
    } else {
        demand.adj <- demand[epoche,]
        natural.adj <- natural.d
    }
    if (sign(natural.adj)== -1) {print("land use cannot be allocated due to spatial restrictions")}
    #
    demand.new <- as.integer(cbind(demand.adj, natural.adj))
#####
#  1.6 Trajectories of land use change
#####
    # general:
    # transitions which are not allowed are set to NA in the respective suitability layer (target)
    # transitions different to 1, referring to transition possible after one iteration (year) are identified 
    # those identified are checked against the transition years vector. if years < transition years the target suitability is set to NA
    # specific steps: 
    # first edit trajectory matrix
    if (length(traj[traj==0 | is.na(traj)] )>0){
      # for not allowed changes (100 years more than modelling years)
      traj[traj==0 | is.na (traj)] <- nrow(demand)+ 100 
    }
#####
#  	1.6.1 Trajectories for land use classes
#####
    # conversion restrictions from all land covers to the  land use classes (suitability layer)
    # for all unique land cover classes to land use classes 

for (i in 1:length(lc_unique)){ 
  #i=5
  # identify classes with restricted trajectories to land use
  traj_ind <-  which(traj[lc_unique[i],lu_suit] != 1) # für 1 an der stelle 2
  #traj_ind
  
  # in case no restriction due to trajectories apply 
  if (length(traj_ind) > 0){  
    # index classes with restricted trajectories
    cat_index <- which(tprop.previous_vector==lc_unique[i])  # für 2 an der stelle 1
    for (a in 1:length(traj_ind)){
      # set p_vector at the specific location for the specific layer  to NA if the amount of years is not reached
      p_vector[ cat_index, traj_ind[a]]<- ifelse (trans.years_vector[cat_index] < traj[lu_suit [traj_ind[a]], lc_unique[i]], NA, p_vector[ cat_index, traj_ind[a]])
    }
  }
}
#tprop.previous_vector
#p_vector


# conversion restrictions from natural vegetation class to any other land cover class
lc.nonatural <-  lc_unique [-(natural)]
for (i in 1:length(lc.nonatural)){
  #i=5
  traj_ind <- which(is.element (1 ,  traj[natural,lc.nonatural[i]])==FALSE) # identify which trajectories are unequal 1 (are not allowed after one year)
  if (length(traj_ind) > 0){ 
    cat_index <- which(tprop.previous_vector==lc.nonatural[i])
    p.natural[cat_index] <- ifelse (trans.years_vector[cat_index] < min(traj[lc.nonatural[i], natural]), NA, p.natural[cat_index])
  }
}
#####

#  1.7 Add elasticities
#####
    # add elas on land use layers
    for (i in 1:length(lu_suit)) {
      elas_index <- which(data_vector==lu_suit[i])
      p_vector[elas_index,i] <- p_vector[elas_index, i] + as.numeric(elas[lu_suit[i]])
    }
    # add ELAS on natural vegetation layer
    for (i in 1: length (natural)){
      elas_index <- which(data_vector== natural[i])
      p.natural[elas_index] <- p.natural[elas_index] + as.numeric(elas[natural[i]])
    }
#####
#  1.8 Combine land use suitability and natural vegetation layer
#####
  p_vector.N <- cbind( p_vector, p.natural)
  #normalize p_vector.N
  for (i in c(1:ncol(p_vector.N))) {
      stretch <- 100/max(p_vector.N[,i],na.rm=TRUE);
      p_vector.N[,i] <- p_vector.N[,i]*stretch;
    }
    
###################################################################################

#####
#  2.3 Run allocation module (defined below)
#####
    if (method == "hierarchical"){
    allocation <- allocation.hierarchical (p_vector.N= p_vector.N,lu.N= lu.N ,demand.new= demand.new, print.log=print.log)
    } 
    if (method == "competitive"){
    allocation <- allocation.module (p_vector.N= p_vector.N,lu.N= lu.N ,demand.new= demand.new, stop.crit= stop.crit, iter.max= iter.max, ncores= ncores, print.plot= print.plot, print.log= print.log)
    }
#####################################################################################
#3. Post-processing
#####
#  3.1 Read results from 2.2
#####
    tprop_vector <- allocation[[1]]
    logfile1  <- rbind (logfile1 , allocation[[2]])    
#####
#  3.2 Spatial restrictions 
#####
#  	3.2.1 Add stable classes (1.4.2)
#####
    tprop_vector[nochange_index] <- data_vector[nochange_index]
#####
#  	3.2.2 Spatial restrictions (1.4.3)    
#####
    tprop_vector[sp.rest_index] <- data_vector[sp.rest_index]
#####
#  3.3 Natural vegetation and succession
#####
#  	3.3.1 Reclassify all natural to pseudo natural classes
#####
    tprop_vector [is.element(tprop_vector, natural)] <- pseudo.N # natural vegetation to pseudo.N class (including areas of spatial restrictions)
##### 
#  	3.3.2 Reclassify pseudo natural class based on trajectory and succession order
#####
    pseudo.index <- which(is.element(tprop_vector, pseudo.N))
	if (length (pseudo.index) >=1){  
    if(length (natural > 1)){
    for (i in 1:length(pseudo.index)){
      #i=1
      #can before.n be translated to natural 
      for (a in length(natural):2){
        if (traj[natural[a], natural[a-1]] < trans.years_vector[pseudo.index[i]] |
              tprop.previous_vector[pseudo.index[i]] == natural[a-1]) {
              tprop_vector[pseudo.index[i]] <- natural[a-1]
        } else{
          tprop_vector[pseudo.index[i]] <- natural[a]
        } 
      }
    }
    } 
	if (length (natural)== 1) {tprop_vector[pseudo.index[i]] <- natural}
	if (sum(is.element (tprop_vector, pseudo.N))!=0) {print( "error in natural vegetation module")}
    }
    # tprop_vector[which(tprop_vector==9)] <- natural[length(natural)]
#####
#  3.4 Saving results and preparing next epoche
#####
#  	3.4.1 Updating transition years vector
#####
    trans.years_vector <- ifelse(tprop_vector==tprop.previous_vector, trans.years_vector + 1, 1) #compare this allocation for transition years, inc if changed, reset to 1 if change
	##write transition years as raster file
	#new.transition <- lc
	#new.transition <- setValues(new.transition, trans.years_vector)
  #assign("global.new.transition", new.transition, envir = .GlobalEnv) 
	#writeRaster(new.transition, paste("transition", epoche, ".tif", sep=""), overwrite=TRUE)
#####
#  	3.4.2 Save final allocation result as raster
#####
    new.data <- lc
    new.data <- setValues(new.data, tprop_vector)
    #assign("global.new.data", new.data, envir = .GlobalEnv) 
    if (print.plot==TRUE){plot(new.data)}
    if (writeRaster==TRUE){writeRaster(new.data, paste("scenario", epoche, ".tif", sep=""), overwrite=TRUE)}
    # name result 
    assign(paste("scenario", epoche, sep=""), new.data)
      
    # now set previous to this epoche and start next
    tprop.previous_vector <- tprop_vector;
    
    print("epoche done")
    #initialize new epoche
    epoche <- epoche+1
  } # end of epoche loop 
######
#4. Return results and logfile
#####
  return(list(stack (mget (paste("scenario", rep(1:nrow(demand)),sep=""))), logfile1))
}


#####################################
# hieracical allocation

allocation.hierarchical <- function (p_vector.N, lu.N , demand.new,  print.log=TRUE) {
  tprop_vector <- rep(NA,times=nrow(p_vector.N) ) 
  demand.init <- demand.new 
  logfile.tmp<- c()
  if  (print.log == TRUE){ pb <- txtProgressBar(min=0, max(length(lu.N)), style=3)}  
  for (i in 1:length (lu.N)) {
    ind <-  order(p_vector.N[,i],decreasing =TRUE ,na.last=NA)# NA's removed
    if (demand.new[order(lu.N)[i]] > length(ind)){
      demand.new[order(lu.N)[i]] <- length(ind)
      print("not enough pixel for allocation")
    }
    tprop_vector [ind [1:demand.new[order(lu.N)[i]]]] <- lu.N[i]
    p_vector.N [ind [1:demand.new[order(lu.N)[i]]],-i] <- NA
    if (print.log == TRUE){
      setTxtProgressBar(pb, i)
    }
  } 
  logfile.tmp <- (demand.new - demand.init)
  if (print.log == TRUE){
    setTxtProgressBar(pb, i)
    close(pb)
    print (logfile.tmp)
  }
  return(list(tprop_vector, logfile.tmp))
}

# allocation <- list(tprop_vector, logfile.tmp)
#writeRaster (new.data, "test_hir4.tif")


#####################################
# competitive allocation
#2. Allocation sub module 
#####
#  2.1 Definition of allocation function
#####
# Description: function providing the allocation routine
#
# p_vector.N | preprocessed suitability vectors including spatial and trajectory based restriction, elasticities etc.
# lu.n | class numbers  of all classes to be modelled (incl. pseudo natural class) 
# demand.new | preprocessed demand file including demand for land use classes plus natural vegetation. Adjusted for spatial restrictions
# stop.crit | stop criteria
# iter.max | maximum of iterations
# ncores | amount of cores
# print.plot | TRUE or FALSE
# print.log | TRUE or FALSE
#####
#  2.2 Function - allocation module
#####
allocation.module <- function(p_vector.N ,lu.N ,demand.new, stop.crit, iter.max, ncores, print.plot, print.log) {
  #####
  #    2.2.1 Initiate & start iteration
  #####
  # initialize iteration
  u = 1
  lu_layer.N <- 1:length(lu.N)# layers of suitability including natural veg. for later use
  lu_suit.N <- length(lu.N) 
  # empty vectors
  logfile1 <- c()
  names.log <- c("u", paste("pix.d",lu_layer.N, sep=""), paste("adj.p",lu_layer.N, sep=""), paste("iter",lu_layer.N, sep=""))
  pix.d.hist <- c()
  perc.d.hist <- c()
  change.p.hist <- c()
  adj.p.hist <- c()
  # iter
  iter <- rep(0,ncol(p_vector.N))
  iter.hist <- c(iter)
  # descriptive variables
  min.demand <- as.integer(which.min(demand.new))
  # initialize cluster
  cl <- makeCluster(getOption("cl.cores", ncores))
  #cl <- makeCluster(ncores, type = "PSOCK") # "MPI" - to test
  on.exit(stopCluster(cl))
  ########
  # start iteration to allocate the requested amount of land use plus natural vegetation
  repeat {
    #####
    #  	2.2.2 Read data & add iter
    ##### 
    if (u==1){
      p2_vector <- p_vector.N
    }else{
      for (i in lu_layer.N){
        p2_vector[,i] <- p_vector.N[,i]+as.numeric(iter[order(lu.N)[i]])
      }
    }
    #####
    #  	2.2.3 Evaluate competitive advantages
    #####
    # main routine to identify competitive advantages between pixels for each location
    
    tprop_vector_tmp <- parRapply(cl,p2_vector,FUN=function(w) ifelse(all(is.na(w)),NA,which.max(w)))
    tprop_vector <- lu.N[tprop_vector_tmp] 
    #####
    #  	2.2.4 Compare amount of received classes to demand     
    #####
    # evaluate result - how many pixels have been assigned to which class
    n <- tabulate(tprop_vector,nbins=max(lu.N))[sort(lu.N)]
    #print(paste("n:" , n))
    # difference of pixels betwen allocated and requested land use/cover
    pix.d  <- n - demand.new 
    perc.d <- pix.d/demand.new
    #save to history
    pix.d.hist <- rbind(pix.d.hist, pix.d)
    perc.d.hist <- rbind(perc.d.hist, perc.d)
    #####      
    #  	2.2.5 Adjust iter for next iteration (starting 2.2.2)
    #####
    if(u==1){ # initializing adj.p 
      adj.p= -1*sign(perc.d) # *1/100
    }else{ # modifying adj.p based on prior results
      change.perc <- abs((pix.d.hist[u-1,]-pix.d.hist[u,])/(pix.d.hist[u-1,])*100) 
      proportion <- abs(pix.d.hist[u,])/ colSums(abs(pix.d.hist[c(u,u-1),]),na.rm=TRUE)
      better<- abs(pix.d.hist[u,])< abs(pix.d.hist[u-1,])
      
      adj.p <- as.vector(ifelse(pix.d.hist[u,]== 0 , 0, # if pixels are correctly allocated no change
                                # if pixels are not correctly allocated 
                                # if last adj hist == 0 then initialize adj.p hist 
                                ifelse(pix.d.hist[u,]!=0 & adj.p.hist[u-1,]==0,-1*sign(perc.d)*1/runif(1,0.5,1.5), # sample(50:150, 1)
                                       #if the last adjustment did not lead to better results and the sign of pix dif is still the same use the double adj.p from last time
                                       ifelse(better==FALSE & sign(pix.d.hist[u,])==sign(pix.d.hist[u-1,]), adj.p.hist[u-1,]*2,
                                              # if the difference is <= 0.001% and <=20 pixel and the sign is the same half of the proportional adjustment              
                                              ifelse(abs(perc.d)<= 0.001 & abs(pix.d.hist[u,])<= 20 & sign(pix.d.hist[u,])==sign(pix.d.hist[u-1,]), (adj.p.hist[u-1,]/2)+((adj.p.hist[u-1,]/2)*proportion),
                                                     # percent change smaller than 20 and the sign is the same, than use the double of last adj.p
                                                     #ifelse(change.perc< 20& sign(perc.d)==sign(perc.d.hist[u-1,]),adj.p.hist[u-1,]*2,
                                                     # percent change is smaller than 40 and larger than 20 - add proportional adj.p to last adj.p
                                                     #ifelse(change.perc< 40& change.perc >= 20 & sign(perc.d)==sign(perc.d.hist[u-1,]), adj.p.hist[u-1,]+(adj.p.hist[u-1,]*proportion),
                                                     ifelse(change.perc< 40&  sign(perc.d)==sign(perc.d.hist[u-1,]), adj.p.hist[u-1,]+(adj.p.hist[u-1,]*proportion),
                                                            # if sign switch adjust proportional  
                                                            ifelse(sign(pix.d.hist[u,])  !=sign(pix.d.hist[u-1,]), -1* adj.p.hist[u-1,]*proportion,
                                                                   # everything else use the last adj.p
                                                                   adj.p.hist[u-1,])))))), mode="numeric") 
      change.p.hist <- rbind(change.p.hist, change.perc)
    }
    #upper and lower boundaries of adj.p
    adj.p <- ifelse (adj.p < -100, -100, ifelse(adj.p > 100,100, adj.p ))
    # adjust iter values for 
    iter <- iter + adj.p
    assign("global.iter", iter , envir = .GlobalEnv) 
    iter <- as.numeric (ifelse(iter < -150, -150, ifelse(iter > 150,150, iter))) # upper and lower bound of iter (should never be reached)
    if (u > 1){
      if (all(sign(iter)==-1) | all(sign(iter)==+1)){ # prevent all iter to have the same sign in the second iteration sign(0) returns 0 
        if (all(sign(iter.hist[nrow(iter.hist),])==-1) | all(sign(iter.hist[nrow(iter.hist),])==+1)){              
          iter[which.min(abs(iter))] <-  0  
          adj.p [which.min(abs(iter))] <- 0 # 
        }}}        
    
    ###
    #save to history
    adj.p.hist <- rbind(adj.p.hist, adj.p)
    iter.hist <- rbind(iter.hist, iter)
    assign("global.iter.hist", iter.hist , envir = .GlobalEnv) 
    iter.hist <-iter.hist
    #####    
    if(print.plot==TRUE){
      plot(0,0,xlim = c(2,iter.max),ylim = c(-100,100),ylab="iter", xlab="iteration", type = "n")
      grid()
      names.legend <- paste ("LC", c (sort(lu.N[-lu_suit.N]),"N"))
      legend("topright", legend=names.legend, col=rainbow(lu_suit.N), pch=15)
      for (i in 1:lu_suit.N){
        lines(c(1:nrow(iter.hist)),iter.hist[,order(lu.N)[i]],col=rainbow(lu_suit.N)[i],type = 'l', lwd=2);
      }
    }
    # write logfile
    log.tmp <- as.vector(c(u, pix.d, adj.p, iter), mode="numeric")
    names(log.tmp) <- names.log
    
    logfile1 <- rbind(logfile1,log.tmp)
    #print tail of logfile
    if(print.log==TRUE){
      #print(tail(logfile1))
      names(log.tmp) <-  c("u", paste ("LC", c (sort(lu.N[-lu_suit.N]),"N")))
      print(log.tmp [1:max(lu_layer.N)+1])
    }
    #####
    #  	2.2.6 Stop iteration when amount of classes ~ demand 
    #####
    #stop argument iteration 
    if (max(abs(perc.d),na.rm=TRUE) < stop.crit[1] & abs(pix.d [min.demand]) <= stop.crit[2] )  {
      break;
    }
    if (abs(pix.d [min.demand]) <= stop.crit[2] & max(abs(pix.d),na.rm=TRUE)<=stop.crit[3]){
      break;
    }
    # stop argument iteration if ITERmax reached and take the ITER with the minimum deviation from the demand from all iterations 
    if (u >= iter.max) {
      current.log <- logfile1[(nrow(logfile1)-iter.max+1):nrow(logfile1),]
      iterfinal_index <-  which.min(rowSums(abs(current.log[,2:(lu_suit.N +1)])))
      iterfinal <- current.log[iterfinal_index, (2*lu_suit.N +2):(3*lu_suit.N +1)]
      if(print.log==TRUE){print(iterfinal)}
      for (i in 1:lu_suit.N){ 
        p2_vector[,i] <- p_vector.N[,i]+ as.numeric(iterfinal[i]);
      }      
      #evaluate best p2_vector
      tprop_vector_tmp <- parRapply(cl,p2_vector,FUN=function(w) ifelse(all(is.na(w)),NA,which.max(w)))
      tprop_vector <- lu.N[tprop_vector_tmp] 
      #
      log.tmp <- as.vector(c(I(iter.max+1), current.log [iterfinal_index,2:(3*lu_suit.N +1)]), mode="numeric")
      logfile1 <- rbind(logfile1, log.tmp)
      if(print.log==TRUE){print("break")
                          print(log.tmp )}
      break;
    }
    #initialize next u sequence
    u=u+1;
  }# next iteration over "u"
  #####
  #  	2.2.7 return allocation vector and logfile
  #####
  #stopCluster(cl);
  print("allocation done")
  return(list(tprop_vector, logfile1))
}




