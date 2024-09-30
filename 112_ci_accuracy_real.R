load("results/forecasts/forecast_true_real.Rdata")
load("dataset_characteristics_real.Rdata")
load("real_datasets.Rdata")

##################################
# CI for individual closed form
##################################

EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}

fc.period <- c(13,26,52)

fc.ind.mean <- fc.ind.q05 <- fc.ind.q25 <- fc.ind.q50 <- fc.ind.q75 <- fc.ind.q95 <- fc.ind.sd <- list()

for (row in 1:22){
  load(paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  purchase_data <- real.datasets[[dataset.characteristics.real$dataset[row]]] 
  
  n <- dataset.characteristics.real$cohort.size[row]
  t <- dataset.characteristics.real$calibration[row]
  
  fc.ind.mean[[row]] <- fc.ind.q05[[row]] <- fc.ind.q25[[row]] <- fc.ind.q50[[row]] <- fc.ind.q75[[row]] <- fc.ind.q95[[row]] <- fc.ind.sd[[row]] <- data.frame("13"=rep(NA,n), "26"=rep(NA,n), "52"=rep(NA,n))
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  
  for (idx in 1:3){
    T.star <- fc.period[idx]
    
    # individual closed form forecast
    x.star.forecast.closed.ind <- matrix(nrow=500,ncol=n)
    for (cust in 1:n)(x.star.forecast.closed.ind[,cust] <- EYt_ind(MCMC_ind[1,,cust],MCMC_ind[2,,cust],MCMC_ind[4,,cust],T.star))
    
    fc.ind.q05[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.05))
    fc.ind.q25[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.25))
    fc.ind.q50[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.5))
    fc.ind.q75[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.75))
    fc.ind.q95[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.95))
    fc.ind.mean[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,mean)
    fc.ind.sd[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,sd)
    
  } # end of fc_rel loop
  rm(individual_data,purchase_data,MCMC_het,params_het, x.star.forecast.closed.het, x.star.forecast.closed.ind,x.star.forecast.mcmc,MCMC_ind)  
  gc()
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop



ci.accuracy.ind.real <- data.frame("50%" = rep(NA,22),"90%" = rep(NA,22), "mean+2sd" = rep(NA,22))
for (row in 1:22){
  ci.accuracy.ind.real[row,1] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.ind.q25[[row]][x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.ind.q75[[row]][x.star.true[[row]][,2] > 0,2]) ) / sum(x.star.true[[row]][,2] > 0)
  ci.accuracy.ind.real[row,2] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.ind.q05[[row]][x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.ind.q95[[row]][x.star.true[[row]][,2] > 0,2]) ) / sum(x.star.true[[row]][,2] > 0)
  ci.accuracy.ind.real[row,3] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= (fc.ind.mean[[row]][x.star.true[[row]][,2] > 0,2]-2*fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2]) & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= (fc.ind.mean[[row]][x.star.true[[row]][,2] > 0,2]+2*fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2])) ) / sum(x.star.true[[row]][,2] > 0)
  
}

# determine relative range to judge on the interval width
ci.range.ind.real <- data.frame("50%" = rep(NA,22),"90%" = rep(NA,22), "mean+2sd" = rep(NA,22))
for (row in 1:22){
  # x.star.true[[row]][,2] > 0
  ci.range.ind.real[row,1] <- mean((fc.ind.q75[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q25[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  ci.range.ind.real[row,2] <- mean((fc.ind.q95[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q05[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  ci.range.ind.real[row,3] <- mean((4* fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  
}

#############################
# CI for simulated mcmc paths
#############################

# Step 1: determine prediction quantiles

source("21_future_transactions_mcmc.R")
load("dataset_characteristics.Rdata")



library(BTYD)

##############################
# FIXED PARAMETERS
##############################


fc.period <- c(13,26,52)




##############################################
# START ALGORITHM
##############################################

# parameter combination loop
#fc.summary <- list()   # already loaded from the first sample
fc.mcmc <- list()

cat("start", format(Sys.time(), '%H:%M:%S'), fill = TRUE)
for (row in 1:22){
  purchase_data <- real.datasets[[dataset.characteristics.real$dataset[row]]] 
  
  load(paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  
  n <- dataset.characteristics.real$cohort.size[row]
  t <- dataset.characteristics.real$calibration[row]
  
  fc.mcmc[[row]] <- array(NA_real_,dim=c(3,n,7))
  dimnames(fc.mcmc[[row]])[[1]] <- fc.period 
  dimnames(fc.mcmc[[row]])[[3]] <- c("q.05", "q.25", "q.50", "q.75","q.95", "mean","sd")
  
  
  
  # create cbs
  cbs <- data.frame("x" = rep(NA,n), "t.x"=rep(NA,n), "T.cal" = rep(NA,n))
  
  
  for (cust in 1:n){
    cbs$x[cust] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs$t.x[cust] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs$T.cal[cust] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  
  for (idx in 1:3){
    T.star <- fc.period[idx]
    
    # MCMC forecast
    x.star.forecast.mcmc <- draw_future_transactions(cbs, MCMC_ind, T.star=T.star, sample_size=500)
    
    
    fc.mcmc[[row]][idx,,1]   <- apply(x.star.forecast.mcmc,2,function(x) quantile(x,0.05))
    fc.mcmc[[row]][idx,,2]   <- apply(x.star.forecast.mcmc,2,function(x) quantile(x,0.25))
    fc.mcmc[[row]][idx,,3]   <- apply(x.star.forecast.mcmc,2,function(x) quantile(x,0.5))
    fc.mcmc[[row]][idx,,4]   <- apply(x.star.forecast.mcmc,2,function(x) quantile(x,0.75))
    fc.mcmc[[row]][idx,,5]   <- apply(x.star.forecast.mcmc,2,function(x) quantile(x,0.95))
    fc.mcmc[[row]][idx,,6]   <- apply(x.star.forecast.mcmc,2, mean)
    fc.mcmc[[row]][idx,,7]   <- apply(x.star.forecast.mcmc,2, sd)
    
    
    save(fc.mcmc,file="results/112a_forecast_mcmc_quantiles_real.Rdata")
  } # end of fc_rel loop
  rm(individual_data,purchase_data,MCMC_het,params_het, x.star.forecast.closed.het, x.star.forecast.closed.ind,x.star.forecast.mcmc,MCMC_ind,params_ind)  
  gc()
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop



# Step 2: determine IC range and accuracy

ci.accuracy.sim.real <- data.frame("50%" = rep(NA,22),"90%" = rep(NA,22), "mean+2sd" = rep(NA,22))
for (row in 1:22){
  ci.accuracy.sim.real[row,1] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,4]) ) / sum(x.star.true[[row]][,2] > 0)
  ci.accuracy.sim.real[row,2] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,1] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,5]) ) / sum(x.star.true[[row]][,2] > 0)
  ci.accuracy.sim.real[row,3] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= (fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,6]-2*fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7]) & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= (fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,6]+2*fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7])) ) / sum(x.star.true[[row]][,2] > 0)
  
}

# determine relative range to judge on the interval width
ci.range.sim.real <- data.frame("50%" = rep(NA,22),"90%" = rep(NA,22), "mean+2sd" = rep(NA,22))
for (row in 1:22){
  # x.star.true[[row]][,2] > 0
  ci.range.sim.real[row,1] <- mean((fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,4] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  ci.range.sim.real[row,2] <- mean((fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,5] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,1]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  ci.range.sim.real[row,3] <- mean((4* fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  
}
save.image("results/112_ci_accuracy_real.Rdata")
