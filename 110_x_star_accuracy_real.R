library(forecast)

load("sample_analysis_real.Rdata")
load("results/forecasts/forecast_MLE_real.Rdata")
load("results/forecasts/forecast_MCMC_real.Rdata")
load("results/forecasts/forecast_heur_real.Rdata")
load("results/forecasts/forecast_true_real.Rdata")

real.idx <- sample.analysis.real$dataset != 1


MAPE.real <- RMSPE.real <- list()
MLE.missing <- vector()

for (T.star in 1:3){
  
  MAPE.real[[T.star]] <- RMSPE.real[[T.star]] <- data.frame("sim_med"=rep(0,22),"ind_med"=rep(0,22),"het_med"=rep(0,22), "sim_mean"=rep(0,22),"ind_mean"=rep(0,22),"het_mean"=rep(0,22),"MLE"=rep(NA,22),"heur"=rep(0,22))
  
  for (row in 1:22){
    
    n <- sample.analysis.real$cohort_size[row]
    t <- sample.analysis.real$calibration[row]
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    RMSPE.real[[T.star]]$sim_med[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,1,2], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star])) #,1,2 = sim, median
    RMSPE.real[[T.star]]$ind_med[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,2,2], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))  #,2,2 = closed ind, median
    RMSPE.real[[T.star]]$het_med[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,3,2], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))  #,3,2 = closed het, median
    
    RMSPE.real[[T.star]]$sim_mean[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,1,1], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))  #,1,1 = sim, mean
    RMSPE.real[[T.star]]$ind_mean[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,2,1], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))  #,2,1 = closed ind, mean 
    RMSPE.real[[T.star]]$het_mean[row] <- 100*(accuracy(x.star.MCMC[[row]][T.star,,3,1], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))  #,3,1 = closed het, mean
    
    if(!is.na(x.star.MLE[[row]][1,T.star])){
    RMSPE.real[[T.star]]$MLE[row] <-  100* (accuracy(x.star.MLE[[row]][,T.star], x.star.true[[row]][,T.star]) [2] / mean(x.star.true[[row]][,T.star]))}  # MLE
    RMSPE.real[[T.star]]$heur[row] <- 100* (accuracy(x.star.heur[[row]][,T.star], x.star.true[[row]][,T.star]) [2] /mean(x.star.true[[row]][,T.star])) # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
   
    MAPE.real[[T.star]]$sim_med[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,1,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,1,2 = sim, median
    MAPE.real[[T.star]]$ind_med[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,2,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,2,2 = closed ind, median
    MAPE.real[[T.star]]$het_med[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,3,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,3,2 = closed het, median
    
    MAPE.real[[T.star]]$sim_mean[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,1,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,1,1 = sim, mean
    MAPE.real[[T.star]]$ind_mean[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,2,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,2,1 = closed ind, mean 
    MAPE.real[[T.star]]$het_mean[row] <- 100*((accuracy(x.star.MCMC[[row]][T.star,,3,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,3,1 = closed het, mean
    
    if(!is.na(x.star.MLE[[row]][1,T.star])){
    MAPE.real[[T.star]]$MLE[row] <- (100*(try((accuracy(x.star.MLE[[row]][,T.star], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star]))))} # MLE
    MAPE.real[[T.star]]$heur[row] <- 100*((accuracy(x.star.heur[[row]][,T.star], x.star.true[[row]][,T.star]) [3])  / mean(x.star.true[[row]][,T.star])) # heur
  
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
   
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}

for (T.star in 1:3){
  RMSPE.real[[T.star]]$MLE <- as.numeric(RMSPE.real[[T.star]]$MLE)
  MAPE.real[[T.star]]$MLE <- as.numeric(MAPE.real[[T.star]]$MLE)
}

save(MAPE.real, RMSPE.real, file="results/110_x_star_deviation_real.Rdata")


