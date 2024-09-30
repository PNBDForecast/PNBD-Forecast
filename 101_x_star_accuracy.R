library(forecast)

load("sample_analysis.Rdata")
load("results/forecasts/forecast_MLE.Rdata")
load("results/forecasts/forecast_MCMC.Rdata")
load("results/forecasts/forecast_heur.Rdata")
load("results/forecasts/forecast_true.Rdata")


MAPE <- RMSPE <- list()
MLE.missing <- vector()

for (T.star in 1:3){
  
  MAPE[[T.star]] <- RMSPE[[T.star]] <- data.frame("sim_med"=rep(0,3000),"ind_med"=rep(0,3000),"het_med"=rep(0,3000), "sim_mean"=rep(0,3000),"ind_mean"=rep(0,3000),"het_mean"=rep(0,3000),"MLE"=rep(NA,3000),"heur"=rep(0,3000))
  
  for (row in 1:3000){
    
    n <- sample.analysis$cohort.size[row]
    t <- sample.analysis$calibration[row]
    
   
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    RMSPE[[T.star]]$sim_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,1,2], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star]) )  #,1,2 = sim, median
    RMSPE[[T.star]]$ind_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,2,2], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star]))  #,2,2 = closed ind, median
    RMSPE[[T.star]]$het_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,3,2], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star])) #,3,2 = closed het, median
    
    RMSPE[[T.star]]$sim_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,1,1], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star])) #,1,1 = sim, mean
    RMSPE[[T.star]]$ind_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,2,1], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star])) #,2,1 = closed ind, mean 
    RMSPE[[T.star]]$het_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,3,1], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star])) #,3,1 = closed het, mean
    
    if (!is.na(x.star.MLE[[row]][1,T.star])){
    RMSPE[[T.star]]$MLE[row] <- 100*((accuracy(x.star.MLE[[row]][,T.star], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star]))   # MLE
    }
    RMSPE[[T.star]]$heur[row] <- 100*((accuracy(x.star.heur[[row]][,T.star], x.star.true[[row]][,T.star]) [2]) / mean(x.star.true[[row]][,T.star])) # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
   
    MAPE[[T.star]]$sim_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,1,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,1,2 = sim, median
    MAPE[[T.star]]$ind_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,2,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,2,2 = closed ind, median
    MAPE[[T.star]]$het_med[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,3,2], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,3,2 = closed het, median
    
    MAPE[[T.star]]$sim_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,1,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,1,1 = sim, mean
    MAPE[[T.star]]$ind_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,2,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,2,1 = closed ind, mean 
    MAPE[[T.star]]$het_mean[row] <- 100*((accuracy(x.star.mcmc[[row]][T.star,,3,1], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) #,3,1 = closed het, mean
    
    if (!is.na(x.star.MLE[[row]][1,T.star])){
    MAPE[[T.star]]$MLE[row] <- 100*((accuracy(x.star.MLE[[row]][,T.star], x.star.true[[row]][,T.star]) [3]) / mean(x.star.true[[row]][,T.star])) # MLE
    }
    MAPE[[T.star]]$heur[row] <- 100*((accuracy(x.star.heur[[row]][,T.star], x.star.true[[row]][,T.star]) [3])  / mean(x.star.true[[row]][,T.star])) # heur
  
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
   
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}

for (T.star in 1:3){
  RMSPE[[T.star]]$MLE <- as.numeric(RMSPE[[T.star]]$MLE)
  MAPE[[T.star]]$MLE <- as.numeric(MAPE[[T.star]]$MLE)
}

save(MAPE, RMSPE, file="results/101_x_star_deviation.Rdata")


boxplot(MAPE[[1]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 13)")
boxplot(MAPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 26)")
boxplot(MAPE[[3]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 52)")

boxplot(MAPE[[1]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 13)")
boxplot(MAPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 26)")
boxplot(MAPE[[3]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 52)")


boxplot(RMSPE[[1]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 13)")
boxplot(RMSPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 26)")
boxplot(RMSPE[[3]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),],outline=F, main="MAPE of future purchases (T* = 52)")

boxplot(RMSPE[[1]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 13)")
boxplot(RMSPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 26)")
boxplot(RMSPE[[3]][unlist(lapply(x.star.true, function(x) sum(x[,3] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 52)")

