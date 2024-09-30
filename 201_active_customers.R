load("results/101_x_star_deviation.Rdata")

##############################################
# define sensitivity and precision functions
##############################################

sensitivity_repeat <- function(x.star, x.star.fc){
  return(100*sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) / sum(x.star > 0,na.rm=TRUE)) # Bezugsgröße: tatsächliche Käufer
}

precision_repeat <- function(x.star, x.star.fc){
  return(100*sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) / sum(x.star.fc > 0.5,na.rm=TRUE)) # Bezugsgröße: als Käufer eingestufte
}

tot.correct.fun <- function(x.star, x.star.fc){
  return(100*((sum(x.star > 0 & x.star.fc >= 0.5,na.rm=TRUE) + sum(x.star == 0 & x.star.fc < 0.5,na.rm=TRUE))/length(x.star))) # Bezugsgröße: als Käufer eingestufte
}


#####################################################
# determine sensitivity, precision, and tot.correct
#####################################################

load("sample_analysis.Rdata")
load("results/forecasts/forecast_MLE.Rdata")
load("results/forecasts/forecast_MCMC.Rdata")
load("results/forecasts/forecast_heur.Rdata")
load("results/forecasts/forecast_true.Rdata")


sensitivity <- precision <-  tot.correct <-  list()

for (T.star in 1:3){
  
  sensitivity[[T.star]] <- precision[[T.star]] <- tot.correct[[T.star]] <- data.frame("sim_med"=rep(NA,3000),"ind_med"=rep(0,3000),"het_med"=rep(0,3000), "sim_mean"=rep(0,3000),"ind_mean"=rep(0,3000),"het_mean"=rep(0,3000),"MLE"=rep(NA,3000),"heur"=rep(0,3000))
  
  for (row in 1:3000){
    
    n <- sample.analysis$cohort.size[row]
    t <- sample.analysis$calibration[row]
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    sensitivity[[T.star]]$sim_med[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,2])   #,1,2 = sim, median
    sensitivity[[T.star]]$ind_med[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    sensitivity[[T.star]]$het_med[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    sensitivity[[T.star]]$sim_mean[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,1])   #,1,1 = sim, mean
    sensitivity[[T.star]]$ind_mean[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    sensitivity[[T.star]]$het_mean[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    sensitivity[[T.star]]$MLE[row] <- as.numeric(try(sensitivity_repeat(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    sensitivity[[T.star]]$heur[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    precision[[T.star]]$sim_med[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,2])   #,1,2 = sim, median
    precision[[T.star]]$ind_med[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    precision[[T.star]]$het_med[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    precision[[T.star]]$sim_mean[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,1])   #,1,1 = sim, mean
    precision[[T.star]]$ind_mean[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    precision[[T.star]]$het_mean[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    precision[[T.star]]$MLE[row] <- as.numeric(try(precision_repeat(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    precision[[T.star]]$heur[row] <- precision_repeat(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    tot.correct[[T.star]]$sim_med[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,2])   #,1,2 = sim, median
    tot.correct[[T.star]]$ind_med[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,2])   #,2,2 = closed ind, median
    tot.correct[[T.star]]$het_med[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,2])   #,3,2 = closed het, median
    
    tot.correct[[T.star]]$sim_mean[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,1,1])   #,1,1 = sim, mean
    tot.correct[[T.star]]$ind_mean[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,2,1])   #,2,1 = closed ind, mean 
    tot.correct[[T.star]]$het_mean[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.mcmc[[row]][T.star,,3,1])   #,3,1 = closed het, mean
    
    tot.correct[[T.star]]$MLE[row] <- as.numeric(try(tot.correct.fun(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star])))  # MLE
    tot.correct[[T.star]]$heur[row] <- tot.correct.fun(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star])  # heur
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}
save(precision, sensitivity, tot.correct, file="results/201_active_customers.Rdata")


##############################################
# plot densities of tot.accuracy
##############################################

T.star.idx <- 3

plot(density(tot.correct[[T.star.idx]]$sim_med), main=paste("tot.correct for T* =",T.star.idx),xlab="",ylab="")
lines(density(tot.correct[[T.star.idx]]$ind_med), col="red")
lines(density(tot.correct[[T.star.idx]]$sim_mean), col="blue")
lines(density(tot.correct[[T.star.idx]]$ind_mean), col="darkgreen")
lines(density(tot.correct[[T.star.idx]]$MLE), col="purple")
lines(density(tot.correct[[T.star.idx]]$het_med), col="gray")
lines(density(tot.correct[[T.star.idx]]$heur), col="orange")

library(psych)
describe(tot.correct[[T.star.idx]])


#####################################################
# determine optimal p threshold and accuracy values
#####################################################


x.star.ind <- lapply(x.star.mcmc, function(x) t(x[,,2,2]))
x.star.sim <- lapply(x.star.mcmc, function(x) t(x[,,1,2]))
rm(x.star.mcmc)

opt.correct.fun <- function(x.star, x.star.fc, p){
  return(100*((sum(x.star > 0 & x.star.fc >= p,na.rm=TRUE) + sum(x.star == 0 & x.star.fc < p,na.rm=TRUE))/length(x.star)))
}


opt.correct <-  opt.p <- list()


p.vec <- seq(0.05,0.95,0.05)

for (T.star in 1:3){
  opt.p[[T.star]] <- data.frame("sim.med"=rep(0,3000),"ind.med"=rep(0,3000),"MLE"=rep(NA,3000),"heur"=rep(0,3000))
  opt.correct[[T.star]] <- data.frame("sim.med"=rep(0,3000),"ind.med"=rep(0,3000),"MLE"=rep(NA,3000),"heur"=rep(0,3000))
  
  for (row in 1:3000){
    
    n <- sample.analysis$cohort.size[row]
    t <- sample.analysis$calibration[row]
    
    temp <- rep(0,19)
    
    # sim.med
    for (p in 1:19){
      temp[p] <- opt.correct.fun(x.star.true[[row]][,T.star], x.star.sim[[row]][,T.star], p.vec[p])   #,1,2 = sim, median
    }
    opt.correct[[T.star]]$sim.med[row] <- max(temp)
    opt.p[[T.star]]$sim.med[row] <- p.vec[which.max(temp)]
    
    # ind.med
    for (p in 1:19){
      temp[p] <- opt.correct.fun(x.star.true[[row]][,T.star], x.star.ind[[row]][,T.star], p.vec[p])   #,1,2 = sim, median
    }
    opt.correct[[T.star]]$ind.med[row] <- max(temp)
    opt.p[[T.star]]$ind.med[row] <- p.vec[which.max(temp)]
    
    # MLE
    for (p in 1:19){
      temp[p] <- as.numeric(try(opt.correct.fun(x.star.true[[row]][,T.star], x.star.MLE[[row]][,T.star], p.vec[p])))  #MLE
    }
    opt.correct[[T.star]]$MLE[row] <- max(temp)
    opt.p[[T.star]]$MLE[row] <- p.vec[which.max(temp)]
    
    # heur
    for (p in 1:19){
      temp[p] <- opt.correct.fun(x.star.true[[row]][,T.star], x.star.heur[[row]][,T.star], p.vec[p])   #heur
    }
    opt.correct[[T.star]]$heur[row] <- max(temp)
    opt.p[[T.star]]$heur[row] <- p.vec[which.max(temp)]
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}
save(opt.p, opt.correct, file="results/201_active_customers_optimised_p.Rdata")






