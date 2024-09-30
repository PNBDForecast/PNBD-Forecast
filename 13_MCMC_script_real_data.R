source("00_general data.R")
load("sample_analysis_real.Rdata")
load("real_datasets.Rdata")



#------------------------

for (row in 1:22){
  
  # 1. Housekeeping I
  MCMC_het <- array(NA_real_,dim=c(500,5))
  dimnames(MCMC_het)[[2]] <- c("r","alpha","s","beta","LL")
  
  purchase_data <- real.datasets[[sample.analysis.real$dataset[row]]]
  
  t <- sample.analysis.real$calibration[row]
  n <- sample.analysis.real$cohort_size[row]
  
  MCMC_ind <- array(NA_real_,dim=c(4,500,n))
  dimnames(MCMC_ind)[[1]] <- c("lambda","mu","tau","z")
  
  cbs <- array(NA_real_,dim=c(n,3))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal")
  
  
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  cbs           <- cbs[1:n,]
  
  
  #  5. define hyper parameters
  hyper_prior <- c("r1"= CV_value^(-2),"r2"= CV_value^(-2)/2.22,
                   "a1"= CV_value^(-2),"a2"= CV_value^(-2)/12.22,
                   "s1"= CV_value^(-2),"s2"= CV_value^(-2)/1,
                   "b1"= CV_value^(-2),"b2"= CV_value^(-2)/17.68)
  
  #  6 perform MCMC
  
  source("10_MCMC_Abe.R")
  helper_array <- array(NA_real_,dim=c(500,n,4))
  for (chain_id in 1:4){
    MCMC_het[((chain_id -1)*mcmc/thin+1):(chain_id*mcmc/thin),] <- alldraws$level_2[[chain_id]]
    for (cust_id in 1:n) {
      helper_array[((chain_id -1)*mcmc/thin+1):(chain_id*mcmc/thin),cust_id,] <- alldraws$level_1[[cust_id]][[chain_id]]
      for (param in 1:4) {
        MCMC_ind[param,,cust_id] <- helper_array[,cust_id,param]
        
      } # end of param loop
    } # end of customer loop 
  }# end of chain loop
  cat("data set",row,format(Sys.time(), '%H:%M:%S'), fill = TRUE)
  
  
  
  save(MCMC_ind, file=paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep=""))
  save(MCMC_het, file=paste("real_est/MCMC_het/het_est_",row,".Rdata",sep=""))
  
  
} # end of row loop
