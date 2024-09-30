
##############################
# FIXED PARAMETERS
##############################

fc.period <- c(13,26,52)

##############################
# HOUSEKEEPING
##############################

load("sample_analysis.Rdata")
x.star.true <- list()


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
  # determine true x.star
  ##########################################
  
  x.star.true[[row]] <- matrix(nrow=n, ncol=3)
  colnames(x.star.true[[row]]) <- fc.period
  
  for (cust in 1:n){
    for (T.star in 1:3){
      x.star.true[[row]][cust,T.star] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*(t+fc.period[T.star])]))-length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))
      }
    }

  cat(format(Sys.time(), '%H:%M:%S'),"Datensatz",row,"beendet.","\n")
  
} # end of row

save(x.star.true,file="results/forecast_true.Rdata")
