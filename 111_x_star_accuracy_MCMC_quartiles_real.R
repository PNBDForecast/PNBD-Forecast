load("results/forecasts/forecast_true_real.Rdata")
load("sample_analysis_real.Rdata")
load("real_datasets.Rdata")

##############################################
# DEFINE INDIVIDUAL FORECAST IN CLOSED FORM
##############################################
EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}

fc.period <- c(13,26,52)

fc.ind.med.q05.real <- fc.ind.med.q95.real <- list()

for (row in 1:22){
  
  purchase_data <- real.datasets[[sample.analysis.real$dataset[row]]]
  load(paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("real_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- sample.analysis.real$cohort_size[row]
  t <- sample.analysis.real$calibration[row]
  
  fc.ind.med.q05.real[[row]] <- fc.ind.med.q95.real[[row]] <- data.frame("13"=rep(NA,n), "26"=rep(NA,n), "52"=rep(NA,n))
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,3))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal")
  
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
    
    fc.ind.med.q05.real[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.05))
    fc.ind.med.q95.real[[row]][,idx] <- apply(x.star.forecast.closed.ind,2,function(x) quantile(x,0.95))
    
  } # end of fc_rel loop
  rm(purchase_data,MCMC_het,params_het, x.star.forecast.closed.het, x.star.forecast.closed.ind,x.star.forecast.mcmc,MCMC_ind,params_ind)  
  gc()
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop


correct.interquartile.ratio.real <- data.frame("13" = rep(NA,22),"26" = rep(NA,22),"52" = rep(NA,22))
for (row in 1:22){
  correct.interquartile.ratio.real[row,1] <- 100*(sum(x.star.true[[row]][,1] >= fc.ind.med.q05.real[[row]][,1] & x.star.true[[row]][,1] <= fc.ind.med.q95.real[[row]][,1]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio.real[row,2] <- 100*(sum(x.star.true[[row]][,2] >= fc.ind.med.q05.real[[row]][,2] & x.star.true[[row]][,2] <= fc.ind.med.q95.real[[row]][,2]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio.real[row,3] <- 100*(sum(x.star.true[[row]][,3] >= fc.ind.med.q05.real[[row]][,3] & x.star.true[[row]][,3] <= fc.ind.med.q95.real[[row]][,3]) ) / nrow(x.star.true[[row]])
}

correct.interquartile.ratio.active.only.real <- data.frame("13" = rep(NA,22),"26" = rep(NA,22),"52" = rep(NA,22))
for (row in 1:22){
  correct.interquartile.ratio.active.only.real[row,1] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,1] > 0,1] >= fc.ind.med.q05.real[[row]][x.star.true[[row]][,1] > 0,1] & x.star.true[[row]][x.star.true[[row]][,1] > 0,1] <= fc.ind.med.q95.real[[row]][x.star.true[[row]][,1] > 0,1]) ) / sum(x.star.true[[row]][,1] > 0)
  correct.interquartile.ratio.active.only.real[row,2] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.ind.med.q05.real[[row]][x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.ind.med.q95.real[[row]][x.star.true[[row]][,2] > 0,2]) ) / sum(x.star.true[[row]][,2] > 0)
  correct.interquartile.ratio.active.only.real[row,3] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,3] > 0,3] >= fc.ind.med.q05.real[[row]][x.star.true[[row]][,3] > 0,3] & x.star.true[[row]][x.star.true[[row]][,3] > 0,3] <= fc.ind.med.q95.real[[row]][x.star.true[[row]][,3] > 0,3]) ) / sum(x.star.true[[row]][,3] > 0)
}

save.image("results/111_90_percent_accuracy_real.Rdata")
