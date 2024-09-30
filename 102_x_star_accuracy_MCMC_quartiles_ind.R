load("results/forecasts/forecast_true.Rdata")
load("sample_analysis.Rdata")

##############################################
# DEFINE INDIVIDUAL FORECAST IN CLOSED FORM
##############################################
EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}

fc.period <- c(13,26,52)

fc.ind.mean <- fc.ind.q05 <- fc.ind.q25 <- fc.ind.q50 <- fc.ind.q75 <- fc.ind.q95 <- fc.ind.sd <- list()

for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
 
  n <- sample.analysis$cohort.size[row]
  t <- sample.analysis$calibration[row]
  
  fc.ind.mean[[row]] <- fc.ind.q05[[row]] <- fc.ind.q25[[row]] <- fc.ind.q50[[row]] <- fc.ind.q75[[row]] <- fc.ind.q95[[row]] <- fc.ind.sd[[row]] <- data.frame("13"=rep(NA,n), "26"=rep(NA,n), "52"=rep(NA,n))
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
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



correct.interquartile.ratio <- data.frame("50%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  correct.interquartile.ratio[row,1] <- 100*(sum(x.star.true[[row]][,2] >= fc.ind.q25[[row]][,2] & x.star.true[[row]][,2] <= fc.ind.q75[[row]][,2]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio[row,2] <- 100*(sum(x.star.true[[row]][,2] >= fc.ind.q05[[row]][,2] & x.star.true[[row]][,2] <= fc.ind.q95[[row]][,2]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio[row,3] <- 100*(sum(x.star.true[[row]][,2] >= (fc.ind.mean[[row]][,2]-2*fc.ind.sd[[row]][,2]) & x.star.true[[row]][,2] <= (fc.ind.mean[[row]][,2]+2*fc.ind.sd[[row]][,2])) ) / nrow(x.star.true[[row]])
  
}

correct.interquartile.ratio.active.only <- data.frame("75%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  correct.interquartile.ratio.active.only[row,1] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.ind.q25[[row]][x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.ind.q75[[row]][x.star.true[[row]][,2] > 0,2]) ) / sum(x.star.true[[row]][,2] > 0)
  correct.interquartile.ratio.active.only[row,2] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.ind.q05[[row]][x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.ind.q95[[row]][x.star.true[[row]][,2] > 0,2]) ) / sum(x.star.true[[row]][,2] > 0)
  correct.interquartile.ratio.active.only[row,3] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= (fc.ind.mean[[row]][x.star.true[[row]][,2] > 0,2]-2*fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2]) & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= (fc.ind.mean[[row]][x.star.true[[row]][,2] > 0,2]+2*fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2])) ) / sum(x.star.true[[row]][,2] > 0)
  
}

# determine relative range to judge on the interval width
relative.range.active <- data.frame("50%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  # x.star.true[[row]][,2] > 0
  relative.range.active[row,1] <- mean((fc.ind.q75[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q25[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.range.active[row,2] <- mean((fc.ind.q95[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q05[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.range.active[row,3] <- mean((4* fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  
}


relative.active.range.75 <- 0
relative.active.range.90 <- 0
relative.active.range.sd <- 0

for (row in 1:3000){
  # x.star.true[[row]][,2] > 0
  relative.active.range.75 <- c(relative.active.range.75,(fc.ind.q75[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q25[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.active.range.90 <- c(relative.active.range.90,(fc.ind.q95[[row]][x.star.true[[row]][,2] > 0,2] - fc.ind.q05[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.active.range.sd <- c(relative.active.range.sd,(4* fc.ind.sd[[row]][x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
}
relative.active.range.75 <- relative.active.range.75[-1]
relative.active.range.90 <- relative.active.range.90[-1]
relative.active.range.sd <- relative.active.range.sd[-1]

save.image("results/102_range_accuracy_ind.Rdata")
