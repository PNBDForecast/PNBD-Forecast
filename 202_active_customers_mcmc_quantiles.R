
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
load("results/22_forecast_mcmc_quantiles.Rdata")
load("results/forecasts/forecast_true.Rdata")


sensitivity <- precision <-  tot.correct <-  list()

for (T.star in 2:2){
  
  sensitivity[[T.star]] <- precision[[T.star]] <- tot.correct[[T.star]] <- data.frame("q.50"=rep(NA,3000),"q.75"=rep(0,3000),"q.95"=rep(0,3000))
  
  for (row in 1:3000){
    
    n <- sample.analysis$cohort.size[row]
    t <- sample.analysis$calibration[row]
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    sensitivity[[T.star]]$q.50[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    sensitivity[[T.star]]$q.75[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    sensitivity[[T.star]]$q.95[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])   
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    precision[[T.star]]$q.50[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    precision[[T.star]]$q.75[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    precision[[T.star]]$q.95[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])   
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    tot.correct[[T.star]]$q.50[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    tot.correct[[T.star]]$q.75[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    tot.correct[[T.star]]$q.95[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}
save(precision, sensitivity, tot.correct, file="results/203_active_customers_mcmc_quantiles.Rdata")




############################################################
# real data sets
############################################################
load("sample_analysis_real.Rdata")
load("results/112_ci_accuracy_real.Rdata")
load("results/forecasts/forecast_true_real.Rdata")

sensitivity.real <- precision.real <-  tot.correct.real <-  list()

for (T.star in 2:2){
  
  sensitivity.real[[T.star]] <- precision.real[[T.star]] <- tot.correct.real[[T.star]] <- data.frame("q.50"=rep(NA,22),"q.75"=rep(0,22),"q.95"=rep(0,22))
  
  for (row in 1:22){
    
    n <- sample.analysis.real$cohort.size[row]
    t <- sample.analysis.real$calibration[row]
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    sensitivity.real[[T.star]]$q.50[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    sensitivity.real[[T.star]]$q.75[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    sensitivity.real[[T.star]]$q.95[row] <- sensitivity_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])   
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    precision.real[[T.star]]$q.50[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    precision.real[[T.star]]$q.75[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    precision.real[[T.star]]$q.95[row] <- precision_repeat(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])   
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    tot.correct.real[[T.star]]$q.50[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,3])   
    tot.correct.real[[T.star]]$q.75[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,4])   
    tot.correct.real[[T.star]]$q.95[row] <- tot.correct.fun(x.star.true[[row]][,T.star], fc.mcmc[[row]][T.star,,5])
    
    #-----------------------------------------------------------------------------------------------------------------------------------------------------
    
    cat(format(Sys.time(), '%H:%M:%S'),"forecast",T.star,"Datensatz",row,"beendet.","\n")
    
  }
  
}

##############################################
# boxplots
##############################################
par(mfrow=c(1,3))
boxplot(tot.correct[[2]], outline=F, main="correctly assigned")
points(rep(1,13), tot.correct.real[[2]][sample.analysis.real$dataset!=1,1],pch=4,col="red")
points(rep(2,13), tot.correct.real[[2]][sample.analysis.real$dataset!=1,2],pch=4,col="red")
points(rep(3,13), tot.correct.real[[2]][sample.analysis.real$dataset!=1,3],pch=4,col="red")

boxplot(sensitivity[[2]], outline=F, main="sensitivity")
points(rep(1,13), sensitivity.real[[2]][sample.analysis.real$dataset!=1,1],pch=4,col="red")
points(rep(2,13), sensitivity.real[[2]][sample.analysis.real$dataset!=1,2],pch=4,col="red")
points(rep(3,13), sensitivity.real[[2]][sample.analysis.real$dataset!=1,3],pch=4,col="red")

boxplot(precision[[2]], outline=F, main="precision")
points(rep(1,13), precision.real[[2]][sample.analysis.real$dataset!=1,1],pch=4,col="red")
points(rep(2,13), precision.real[[2]][sample.analysis.real$dataset!=1,2],pch=4,col="red")
points(rep(3,13), precision.real[[2]][sample.analysis.real$dataset!=1,3],pch=4,col="red")


