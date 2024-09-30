#####################################
# sensitivity.real and precision.real functions
#####################################


sensitivity.real_repeat <- function(x.star, x.star.fc){
  return(100*(sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) / sum(x.star > 0,na.rm=TRUE))) # Bezugsgröße: tatsächliche Käufer
}

precision.real_repeat <- function(x.star, x.star.fc){
  return(100*(sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) / sum(x.star.fc > 0.5,na.rm=TRUE))) # Bezugsgröße: als Käufer eingestufte
}

tot.correct.real.fun <- function(x.star, x.star.fc){
  return(100*((sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) + sum(x.star == 0 & x.star.fc < 0.5,na.rm=TRUE))/length(x.star))) # Bezugsgröße: als Käufer eingestufte
}


load("sample_analysis_real.Rdata")
load("results/forecasts/forecast_MLE_real.Rdata")
load("results/forecasts/forecast_MCMC_real.Rdata")
load("results/forecasts/forecast_heur_real.Rdata")
load("results/forecasts/forecast_true_real.Rdata")


sensitivity.real <- precision.real <-  tot.correct.real <- list()

for (T.star in 1:3){
  
  sensitivity.real[[T.star]] <- tot.correct.real[[T.star]] <- precision.real[[T.star]] <- data.frame("sim_med"=rep(0,22),"ind_med"=rep(0,22),"het_med"=rep(0,22), "sim_mean"=rep(0,22),"ind_mean"=rep(0,22),"het_mean"=rep(0,22),"MLE"=rep(NA,22),"heur"=rep(0,22))
  
  for (row in 1:22){
    
    n <- sample.analysis.real$cohort.size[row]
    t <- sample.analysis.real$calibration[row]
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    sensitivity.real[[T.star]]$sim_med[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,2])   #,1,2 = sim, median
    sensitivity.real[[T.star]]$ind_med[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    sensitivity.real[[T.star]]$het_med[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    sensitivity.real[[T.star]]$sim_mean[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,1])   #,1,1 = sim, mean
    sensitivity.real[[T.star]]$ind_mean[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    sensitivity.real[[T.star]]$het_mean[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    sensitivity.real[[T.star]]$MLE[row] <- as.numeric(try(sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    sensitivity.real[[T.star]]$heur[row] <- sensitivity.real_repeat(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    precision.real[[T.star]]$sim_med[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,2])   #,1,2 = sim, median
    precision.real[[T.star]]$ind_med[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    precision.real[[T.star]]$het_med[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    precision.real[[T.star]]$sim_mean[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,1])   #,1,1 = sim, mean
    precision.real[[T.star]]$ind_mean[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    precision.real[[T.star]]$het_mean[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    precision.real[[T.star]]$MLE[row] <- as.numeric(try(precision.real_repeat(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    precision.real[[T.star]]$heur[row] <- precision.real_repeat(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    tot.correct.real[[T.star]]$sim_med[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,2])   #,1,2 = sim, median
    tot.correct.real[[T.star]]$ind_med[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    tot.correct.real[[T.star]]$het_med[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    tot.correct.real[[T.star]]$sim_mean[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,1,1])   #,1,1 = sim, mean
    tot.correct.real[[T.star]]$ind_mean[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    tot.correct.real[[T.star]]$het_mean[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MCMC[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    tot.correct.real[[T.star]]$MLE[row] <- as.numeric(try(tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    tot.correct.real[[T.star]]$heur[row] <- tot.correct.real.fun(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}

save(precision.real,sensitivity.real, tot.correct.real , file="results/212_active_customers_real.Rdata")



#####################################################
# determine optimal p threshold and accuracy values
#####################################################


x.star.ind <- lapply(x.star.MCMC, function(x) t(x[,,2,2]))
x.star.sim <- lapply(x.star.MCMC, function(x) t(x[,,1,2]))
rm(x.star.MCMC)

opt.correct.real.fun <- function(x.star, x.star.fc, p){
  return(100*((sum(x.star > 0 & x.star.fc >= p,na.rm=TRUE) + sum(x.star == 0 & x.star.fc < p,na.rm=TRUE))/length(x.star)))
}


opt.correct.real <-  opt.p.real <- list()


p.vec <- seq(0.05,0.95,0.05)

for (T.star in 1:3){
  opt.p.real[[T.star]] <- data.frame("sim.med"=rep(0,22),"ind.med"=rep(0,22),"MLE"=rep(NA,22),"heur"=rep(0,22))
  opt.correct.real[[T.star]] <- data.frame("sim.med"=rep(0,22),"ind.med"=rep(0,22),"MLE"=rep(NA,22),"heur"=rep(0,22))
  
  for (row in 1:22){
    
    n <- sample.analysis.real$cohort.size[row]
    t <- sample.analysis.real$calibration[row]
    
    temp <- rep(0,19)
    
    # sim.med
    for (p in 1:19){
      temp[p] <- opt.correct.real.fun(x.star.true[[row]][,T.star], x.star.sim[[row]][,T.star], p.vec[p])   #,1,2 = sim, median
    }
    opt.correct.real[[T.star]]$sim.med[row] <- max(temp)
    opt.p.real[[T.star]]$sim.med[row] <- p.vec[which.max(temp)]
    
    # ind.med
    for (p in 1:19){
      temp[p] <- opt.correct.real.fun(x.star.true[[row]][,T.star], x.star.ind[[row]][,T.star], p.vec[p])   #,1,2 = sim, median
    }
    opt.correct.real[[T.star]]$ind.med[row] <- max(temp)
    opt.p.real[[T.star]]$ind.med[row] <- p.vec[which.max(temp)]
    
    # MLE
    for (p in 1:19){
      temp[p] <- as.numeric(try(opt.correct.real.fun(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star], p.vec[p])))  #MLE
    }
    opt.correct.real[[T.star]]$MLE[row] <- max(temp)
    opt.p.real[[T.star]]$MLE[row] <- p.vec[which.max(temp)]
    
    # heur
    for (p in 1:19){
      temp[p] <- opt.correct.real.fun(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star], p.vec[p])   #heur
    }
    opt.correct.real[[T.star]]$heur[row] <- max(temp)
    opt.p.real[[T.star]]$heur[row] <- p.vec[which.max(temp)]
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}
save(opt.p.real, opt.correct.real, file="results/211_active_customers_optimised_p_real.Rdata")
