source("00_general data.R")
source("21_future_transactions_mcmc.R")
load("dataset_characteristics.Rdata")

fc.summary <- list()


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
for (row in 199:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- dataset.characteristics$cohort.size[row]
  t <- dataset.characteristics$calibration[row]

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
  
  
  for (idx in 3:3){
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
      
      
    
  } # end of fc_rel loop
    rm(individual_data,purchase_data,MCMC_het,params_het, x.star.forecast.closed.het, x.star.forecast.closed.ind,x.star.forecast.mcmc,MCMC_ind,params_ind)  
    gc()
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop

save(fc.mcmc,file="results/22_forecast_mcmc_quantiles_52.Rdata")
