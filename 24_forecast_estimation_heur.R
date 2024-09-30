
##############################
# FIXED PARAMETERS
##############################

fc.period <- c(13,26,52)

##############################
# HOUSEKEEPING
##############################

load("sample_analysis.Rdata")
x.star.heur <- list()


for (row in 1:3000){
  
  ##########################################
  # PREPARE DATA SET
  ##########################################
  
  load(paste("sim_data/dataset_",row,".Rdata",sep=""))
  n <- sample.analysis$cohort.size[row]
  t <- sample.analysis$calibration[row]
  
    # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  
    
  ##########################################
  # heuristic forecast
  ##########################################
  
  x.star.heur[[row]] <- matrix(nrow=n, ncol=3)
  colnames(x.star.heur[[row]]) <- fc.period
  
  for (T.star in 1:3){
    x.star.heur[[row]][,T.star] <- cbs[,1]/t * fc.period[T.star]
    }

  cat(format(Sys.time(), '%H:%M:%S'),"Datensatz",row,"beendet.","\n")
  
} # end of row

save(x.star.heur,file="results/forecast_heur.Rdata")
