
# We only consider the individual closed form, MLE and heur here.
library(BTYD)
source("00_general data.R")
load("sample_analysis.Rdata")
load("sim_est/MLE_est_CLVTools.Rdata")

T.star <- 130

fc.ind <- fc.MLE <- fc.heur <- fc.true <- fc.ind.mean <-  list()


# DEFINE INDIVIDUAL FORECAST IN CLOSED FORM
EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}

# start algorithm

for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- sample.analysis$cohort.size[row]
  t <- sample.analysis$calibration[row]
  
  fc.ind[[row]] <- fc.MLE[[row]] <- fc.heur[[row]] <- fc.true[[row]] <- vector(length=n)
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  # #-------------------------------------
  # individual closed form forecast
  temp <- matrix(nrow=500,ncol=n)
  for (cust in 1:n)(temp[,cust] <- EYt_ind(MCMC_ind[1,,cust],MCMC_ind[2,,cust],MCMC_ind[4,,cust],T.star))
  fc.ind[[row]] <- apply(temp,2,median)
  
  # individual closed mean
  temp <- matrix(nrow=500,ncol=n)
  for (cust in 1:n)(temp[,cust] <- EYt_ind(MCMC_ind[1,,cust],MCMC_ind[2,,cust],MCMC_ind[4,,cust],T.star))
  fc.ind.mean[[row]] <- apply(temp,2,mean)

  #-------------------------------------
  # MLE
  fc.MLE[[row]] <- pnbd.ConditionalExpectedTransactions(MLE.est[row,1:4],T.star=T.star,cbs[,1],cbs[,2],cbs[,3],hardie = TRUE)

  #-------------------------------------
  # heur
  fc.heur[[row]] <- cbs[,1]/cbs[,3] * T.star

  #-------------------------------------
  # true
  for (cust in 1:n){
  fc.true[[row]][cust] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*(t+T.star)]))-length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))
  }

rm(individual_data,purchase_data,MCMC_het)  
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop

save(fc.heur, fc.ind, fc.MLE, fc.true, fc.ind.mean, file="results/forecasts_topA.Rdata")
