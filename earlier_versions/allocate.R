 #
# file:     allocate.R
#
# coder:    moenkemt@geo.hu-berlin.de
#			florian.gollnow@geo.hu-berlin.de
#
# purpose:  iterate V to achieve demand from TC
#           variable names according to VERBURG2002
#


# load libs
library("parallel")#, lib.loc="/usr/lib/R/library")
library("sp")#, lib.loc="/usr/lib/R/library")
library("raster")#, lib.loc="/usr/lib/R/library")
library("rgdal")#, lib.loc="/usr/lib/R/library")


# load list of files as rasters and return as matrix
load_rasters_as_matrix = function (filenames) {
  P=dim(0)
  for (rasterfile in filenames) {
    thisfile = as.character(rasterfile)
    print(thisfile)
    if (is.na(thisfile)) {
      # empty filename, do not load, fill with zero
      P=rbind(P,null_vector);
    } else {
      # filename, load raster and add values to P
      P=rbind(P,getValues(raster(thisfile)));
    }
  }
  return(as.matrix(P));
}

normalize_prob_vectors=function(P) {
  layer=nrow(P)
  for (i in c(1:layer)) {
    stretch <- 1/max(P[i,],na.rm=TRUE);
    P[i,] <- P[i,]*stretch;
  }
  return(P)
}

plot_ITER_values=function(ITERhist,layer) {
  plot(0,0,xlim = c(2,ITERmax-10),ylim = c(-1.0,1.0),ylab="ITER", xlab="iteration", type = "n")
  for (i in c(1:layer)) {
    lines(c(1:nrow(ITERhist)),ITERhist[,i],col = palette[i],type = 'l', lwd=1.5);
  }
  legend("topright", c("forest", "secondary", "pasture", "cropland", "O", "O", "urban", "other"),col=palette, lty=1 , cex=0.8)
}




# knobs and switches
epoche <- 1;
ITERmax <- 210; # maximum of iteration in each epoche 200

ELAS <- c(0.3, 0.2, 0.7, 0.7, 0.0, 0.0, 1.0, 0.7) # 8 ELAS values

# load modell configuration
setwd("~/AlucR")
demand=read.csv2("demand_linear2.csv")
trajectories=read.csv("trajectories.csv", header=FALSE)
palette <- c("darkgreen", "lightgreen", "goldenrod4", "gold", "grey90", "blue","red", "grey")

# prepare logdata
logdata_all <- c("epoche","u","forestperc","secondperc","pastureperc","agriperc","waterperc","maskperc","urbanperc","otherperc","meanperc","forestdiff","seconddiff","pasturediff","agridiff","waterdiff","maskdiff","urbandiff","otherdiff", 
                 "forestIter","secondIter","pastureIter","agriIter","waterIter","maskIter","urbanIter","otherIter")
logdata=c("epoche","u","forestperc","secondperc","pastureperc","agriperc","waterperc","maskperc","urbanperc","otherperc","meanperc","forestdiff","seconddiff","pasturediff","agridiff","waterdiff","maskdiff","urbandiff","otherdiff")


# start loop over epoches in demand CSV
repeat {
  # load initial data
  data_raster <- raster(as.character(demand$data[epoche]))
  data_vector <- getValues(data_raster)

  # prepare P_matrix
  null_vector <- c(rep(0,length(data_vector)));
  P_vector=load_rasters_as_matrix(demand[epoche,3:10]);
  protected_vector <- getValues(raster(as.character(demand$protected[epoche])));
  print("reading main data - done")
  
  # mask classes not included in the model, set to NA all we dont need
  #data_vector[which (data_vector==0 | data_vector==5 | data_vector==6 | data_vector>=8 ) ] <- NA 
  data_vector[which (data_vector==0 | data_vector==9 ) ] <- NA 
  
  # also do this for our rasters in P matrix
  P_vector [, is.na (data_vector) ] <- NA
  print("masking NA - done")
  
  # data raster (TC) & subset only in first epoche
	if (epoche==1) { 
    # for the first epoche the previous TPROP is our data
	  TPROP_previous_vector<-data_vector;
	  transition_years_vector<-rep(5,length(TPROP_previous_vector))# 5 years under the current land use for initiation  - any better ideas?
    
  #prepare index for buffer of urban areas
  #data_urban_raster=data_raster; 
  #now all urban areas to 7, otherwise NA
  #data_urban_raster[data_urban_raster!=7]=NA;
  #data_urban_buffer_raster=buffer(data_urban_raster,width=2000);
  #writeRater(data_urban_buffer_raster, "urban_buffer.tif")
  
  #prepared buffer from prepared 
  data_urban_buffer_raster <- raster("urban_buffer.tif")  
  
  data_urban_buffer_vector <- getValues(data_urban_buffer_raster);
  # which areas are out of urban area buffer?
  data_urban_buffer_index <- which(is.na(data_urban_buffer_vector));
  print("buffer index - done - only first iteration")
  
	}
  

  # copy dimension of P_vector to work with
  P2_vector <- P_vector;   

  # initilize TPROP
  TPROP_vector <- null_vector;

  # get dimensions
  layer=nrow(P_vector)
  len=ncol(P_vector)
  
  # algorithmus:
  # which pixel we need to examine at all?
  TPROP_previous_landuse=which(!is.na(TPROP_previous_vector))
  # which layers are used the previous epoche?
  TPROP_previous_layers=as.vector(na.exclude(unique(TPROP_previous_vector)))
  # schleife i in 1-length(TRPOP_previous_vector)
  for (i in TPROP_previous_landuse) {
    # enthält 1 bis 8 als pixel auf voriger epoche
    # lu=TPROP_prev_vector[i]
    lu=TPROP_previous_vector[i];
    # lu ist dann die landnutzung an dieser stelle i
    # dann: schleife t in 1 bis length(trajectories[lu,])
    for (t in TPROP_previous_layers) {
      # t ist quasi die potentiell neue landnutzung
      # wenn trajectories[lu,t] noch nicht transistion_years
      # setzen wir P_vector(t,i) auf NA
      if (transition_years_vector[i]<trajectories[lu,t]) {
        P_vector[t,i]=NA;
      }
    }
  }
  print("trajectories/matrix - done")  
    
    
  # pixel where was forest, or secondary for >10 years
  #forest_transition_index=which((TPROP_previous_vector==1) | ((TPROP_previous_vector==2) & (transition_years_vector>10)));
  # all other set to NA in P forest
  #P_vector[1,-(forest_transition_index)] <- NA;
  
  # No direct conversion of forest to secondary
  #secondary_transition_index <- which(TPROP_previous_vector==1)
  #P_vector[2,secondary_transition_index] <- NA;
  
  # prepare cropland/agri areas
  #cropland_transition_index <- which(TPROP_previous_vector!=1)
  # no direct conversion of forest to cropland
  #P_vector[4,!cropland_transition_index] <- NA;

  # set P to NA where no urban possible
  #P_vector[7,data_urban_buffer_index] <- NA;
  
  # print("trajectories - done")
  
  # normalize raster vector
  P_vector=normalize_prob_vectors(P_vector)
  print ("stretch - done")
 
  # Spatial restrictions
  protected_index <- which(!is.na(protected_vector));

  for (i in c(1:layer)) {
    P_vector[i,protected_index] <- NA;
    }
  
  # read Demand
  
  D_data <- c(demand$d_forest[epoche],demand$d_secondary[epoche],demand$d_pasture[epoche],demand$d_agri[epoche],demand$d_mask[epoche],demand$d_water[epoche],demand$d_urban[epoche],demand$d_other[epoche])
  D_protected <- tabulate(data_vector[!is.na(protected_vector)],nbins=layer) 
  
  D <- D_data-D_protected; 
  
  # add ELAS auf P for every layer
  for (i in c(1:layer)) {
    ELAS_index <- which(data_vector==i)
    P_vector[i,ELAS_index] <- P_vector[i,ELAS_index]+ELAS[i];
  }
  
  # set initial values for ITER 
  if (epoche==1) {
    ITER <- rep(0,layer);
  } else {
    ITER <- ITERhist[nrow(ITERhist),]; #from the last epoche  
  }
  
	ITERhist <- rbind(ITER)
	#empty objects for history
	DiffPixhist <- NULL # 
	Diffhist <- NULL
	ChangePHist <- NULL
	KorrHist <- NULL 
  	
	#start iteration
	print ("start iteration")
	u=10; # start ITER loop
	#make cluster 
	cl <- makeCluster(getOption("cl.cores", detectCores()/2));
  repeat { 
    # calc P2 values with ITER 
    if (u==10 & epoche==1){
      P2_vector<- P_vector
    }else{
    for (i in c(1:layer)) {
      P2_vector[i,] <- P_vector[i,]+ITER[i];
    }
	}
    
    # find winner
    # if last layer is NA, dann NA, sonst which.max(P2)
    TPROP_vector <- parCapply(cl,P2_vector,FUN=function(w) ifelse(is.na(w[8]),NA,which.max(w)));
       
    # find distribution of this iteration as N
    N <- tabulate(TPROP_vector,nbins=layer);
    
   
    # calculate percentage of wrong allocated pixels per class
    Diff_Pix <- N-D
    Diff <- Diff_Pix/D;
    
    # cut to 100% 
    Diff <- ifelse(Diff>1,1,ifelse (Diff< -1,-1, Diff))
	
	#save history
	DiffPixhist <- rbind (DiffPixhist, Diff_Pix)
    Diffhist <- rbind (Diffhist, Diff)
    

	if(u==10){ # initializing Korr 
	Korr= -1*sign(Diff)*1/100
	}else{ # modifying Korr based on prior results
	ChangeP <- abs((DiffPixhist[u-10,]-DiffPixhist[u-9,])/(DiffPixhist[u-10,])*100) 
	proportion <- abs(DiffPixhist[u-9,])/ colSums(abs(DiffPixhist[c(u-9,u-10),]),na.rm=TRUE)
	better<- abs(DiffPixhist[u-9,])< abs(DiffPixhist[u-10,])
  
	Korr <- ifelse(DiffPixhist[u-9,]		== 0 , 																0,
			ifelse(DiffPixhist[u-9,]		!=0 	& KorrHist[u-10,]			==0,							-1*sign(Diff)*1/100,
            ifelse(better					==FALSE & sign(DiffPixhist[u-9,])	==sign(DiffPixhist[u-10,]), 	KorrHist[u-10,]*2,
            ifelse(abs(Diff) 				<= 0.001 & abs(DiffPixhist[u-9,])<= 20 & sign(DiffPixhist[u-9,])	==sign(DiffPixhist[u-10,]), 	(KorrHist[u-10,]/2)+(KorrHist[u-10,]/2)*proportion,
			ifelse(ChangeP 					< 20 	& sign(Diff)				==sign(Diffhist[u-10,]), 		KorrHist[u-10,]*2,
            ifelse(ChangeP 					< 40 	& ChangeP >= 20 & sign(Diff)==sign(Diffhist[u-10,]), 		KorrHist[u-10,]+(KorrHist[u-10,]*proportion),
            ifelse(sign(DiffPixhist[u-9,])	!=sign(DiffPixhist[u-10,]),  										-1* KorrHist[u-10,]*proportion,
																												KorrHist[u-10,])))))))
	
	#print conditions
	print(paste(u, epoche, sep=":"))
	cat("1:", DiffPixhist[u-9,]== 0)
	print(paste(u, epoche, sep=":"))
	cat("2:",better==FALSE & sign(DiffPixhist[u-9,])==sign(DiffPixhist[u-10,]))
	print(paste(u, epoche, sep=":"))
	cat("3:",ChangeP < 20 & sign(Diff)==sign(Diffhist[u-10,]))  
	print(paste(u, epoche, sep=":"))
	cat("4:", ChangeP < 40 & ChangeP >= 20 & sign(Diff)==sign(Diffhist[u-10,]))  
	print(paste(u, epoche, sep=":"))
	cat("5:", sign(DiffPixhist[u-9,])!=sign(DiffPixhist[u-10,]) )
    
    ChangePHist <- rbind(ChangePHist, ChangeP)
    }
	# and apply to ITER
    ITER <- ITER + Korr;
    
	#save History
	KorrHist <- rbind(KorrHist, Korr)
	ITERhist <- rbind(ITERhist,ITER);
    
	# save ITER for the plot and show ITER values
  plot_ITER_values(ITERhist,layer);

  #print main info
	print(paste(u, epoche, sep=":")) # cat() same to print  
	cat("Diff:",Diff);
	print(paste(u, epoche, sep=":"))
	cat("ITER:",ITER);
	print(paste(u, epoche, sep=":"))
	cat("Korr:",Korr)
	print(paste(u, epoche, sep=":"))
	cat("Pixel_diff:",Diff_Pix)
	print(paste(u, epoche, sep=":"))
  
	print(tail(DiffPixhist))
	print(tail(KorrHist))

	###
    # add to logdata_all  
    logdata_all <- rbind(logdata_all,c(epoche,u,Diff,mean(Diff,na.rm=TRUE),Diff_Pix, ITER));
    
    #initialize next u sequence
    u=u+1;
    
    # stop if max Diff <0.3% & urban Pixel diff= +-1 
    #if (max(abs(Diff),na.rm=TRUE) < 0.001 & abs(Diff_Pix [7]) <= 1 & max(abs(Diff_Pix),na.rm=TRUE)<=10) {
    #  break;
    #}
  
	if (max(abs(Diff),na.rm=TRUE) < 0.001 & abs(Diff_Pix [7]) <= 1 & max(abs(Diff_Pix),na.rm=TRUE)<=10)  {
	    break;
  	}
  if (abs(Diff_Pix [7]) <= 1 & max(abs(Diff_Pix),na.rm=TRUE)<=10){
    break;
  }
	
    # stop iteration if ITERmax reached
    if (u > ITERmax) { 
      break;
    }
  
  } # next iteration over "u"
  #stop cluster
  stopCluster(cl);
  print("allocation done")
  
  # add logdata of final allocation for the epoche 
  logdata <- rbind(logdata,c(epoche,u,Diff,mean(Diff,na.rm=TRUE),Diff_Pix))
 
  # Add PA back from data
  TPROP_vector[protected_index] <- data_vector[protected_index];
  
  # fill water back from data
  data_water_index <- which(data_vector==6);
  TPROP_vector[data_water_index] <- 6;
  
  # save raster for new epoche
  epoche <- epoche+1
  
  winner_raster <- data_raster
  winner_raster <- setValues(winner_raster, TPROP_vector)
  #plot(winner_raster)
  winner_filename <- as.character(demand$data[epoche])
  if (winner_filename=="") {	
    print("all demand done.")
    break;
  }
  # write and plot/save raster
  writeRaster(winner_raster, winner_filename, datatype="INT1U", overwrite=TRUE);
  plot(winner_raster, col=palette);
  
  # write logdata to file
  write.table(logdata_all,file="logdata_all_25_11_GEOS2_v1.csv",col.names=FALSE,sep=",");
  write.table(logdata,file="logdata_25_11_GEOS2_v1.csv",col.names=FALSE,sep=",");
  
  # compare this allocation for transistion years, inc if changed, reset to 1 if change
  transition_years_vector <- ifelse(TPROP_vector==TPROP_previous_vector,transition_years_vector+1,1);
 
  # now set previous to this epoche and start next
  TPROP_previous_vector <- TPROP_vector;

  print("epoche done")
} # end of epoche loop

