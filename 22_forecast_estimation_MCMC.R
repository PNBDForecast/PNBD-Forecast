source("00_general data.R")
source("21_future_transactions_mcmc.R")
load("2023_dataset_sample.Rdata")

fc.summary <- list()


library(BTYD)

##############################
# FIXED PARAMETERS
##############################


fc.period <- c(13,26,52)


##############################################
# DEFINE INDIVIDUAL FORECAST IN CLOSED FORM
##############################################
EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}


##############################################
# START ALGORITHM
##############################################

# parameter combination loop
#fc.summary <- list()   # already loaded from the first sample

cat("start", format(Sys.time(), '%H:%M:%S'), fill = TRUE)
for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- sample.analysis$cohort.size[row]
  t <- sample.analysis$calibration[row]

  fc.summary[[row]] <- array(NA_real_,dim=c(3,n,3,2))
  dimnames(fc.summary[[row]])[[1]] <- fc.period 
  dimnames(fc.summary[[row]])[[3]] <- c("simulation","closed_ind", "closed_het")
  dimnames(fc.summary[[row]])[[4]] <- c("mean","median")
      
  
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

      # MCMC forecast
      x.star.forecast.mcmc <- draw_future_transactions(cbs, MCMC_ind, T.star=T.star, sample_size=500)
      
      # individual closed form forecast
      x.star.forecast.closed.ind <- matrix(nrow=500,ncol=n)
      for (cust in 1:n)(x.star.forecast.closed.ind[,cust] <- EYt_ind(MCMC_ind[1,,cust],MCMC_ind[2,,cust],MCMC_ind[4,,cust],T.star))
      
      # het closed form forecast
      x.star.forecast.closed.het <- matrix(nrow=500,ncol=n)
      for (draw in 1:nrow(x.star.forecast.closed.het)) (x.star.forecast.closed.het[draw,] <- pnbd.ConditionalExpectedTransactions(MCMC_het[draw,1:4],T.star=T.star,cbs[,1],cbs[,2],cbs[,3],hardie = TRUE))
      
      ##########################################
      # SUMMARISE FORECASTS
      ##########################################
      
      ### fc.summary[,,,,1]: MCMC simulation
      fc.summary[[row]][idx,,1,1]   <- apply(x.star.forecast.mcmc,2,mean) 
      fc.summary[[row]][idx,,1,2]   <- apply(x.star.forecast.mcmc,2,median)
      
      ### fc.summary[,,,,2]: individual closed form
      fc.summary[[row]][idx,,2,1]   <- apply(x.star.forecast.closed.ind,2,mean) 
      fc.summary[[row]][idx,,2,2]   <- apply(x.star.forecast.closed.ind,2,median)
     
      ### fc.summary[,,,,3]: het closed form
      fc.summary[[row]][idx,,3,1]   <- apply(x.star.forecast.closed.het,2,mean) 
      fc.summary[[row]][idx,,3,2]   <- apply(x.star.forecast.closed.het,2,median)
      
      
    save(fc.summary,file="results/forecast_estimates.Rdata")
  } # end of fc_rel loop
    rm(individual_data,purchase_data,MCMC_het,params_het, x.star.forecast.closed.het, x.star.forecast.closed.ind,x.star.forecast.mcmc,MCMC_ind,params_ind)  
    gc()
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop

