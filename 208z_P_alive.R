source("00_general data.R")
source("21_future_transactions_mcmc.R")
load("2023_dataset_sample.Rdata")
load("results/forecast_true.Rdata")
load("results/MLE_est_CLVtools.Rdata")
library(CLVTools)

#######################################
# Step 1: determine p(alive)
#######################################


p.alive.mcmc.median <- p.alive.ind.median <- p.alive.mcmc.mean <- p.alive.ind.mean <- p.alive.MLE <- list()

P.alive.ind <- function(t.x, T, lambda, mu){
  return((1+(lambda / (lambda + mu)) * (exp((lambda + mu)*(T - t.x)) - 1))^(-1))
}

for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- dataset.sample$cohort_size[row]
  t <- dataset.sample$calibration[row]
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop

  p.alive.mcmc.median[[row]] <- apply(MCMC_ind[4,,], 2,median)
  p.alive.mcmc.mean[[row]] <- apply(MCMC_ind[4,,], 2,mean)

  p.alive.ind.median[[row]] <- sapply(1:n, function(x) median(P.alive.ind(cbs[x,2], cbs[x,3], MCMC_ind[1,,x], MCMC_ind[2,,x])))
  p.alive.ind.mean[[row]] <- sapply(1:n, function(x) mean(P.alive.ind(cbs[x,2], cbs[x,3], MCMC_ind[1,,x], MCMC_ind[2,,x])))

  p.alive.MLE[[row]] <- CLVTools:::pnbd_nocov_PAlive(MLE.est[row,1],MLE.est[row,2],MLE.est[row,3],MLE.est[row,4], cbs[,1], cbs[,2], cbs[,3])
  
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop

save(p.alive.ind.mean, p.alive.ind.median, p.alive.mcmc.mean, p.alive.mcmc.median, p.alive.MLE, file="results/P_alives.Rdata")

################################################
### step 2: derive activity status from p(alive)
#################################################

p.vec <- seq(0,1,0.05)

activity.mcmc.median <- activity.mcmc.mean <- activity.ind.median <- activity.ind.mean <- activity.MLE <- data.frame(rep(NA, 3000))

for (i in 1:21){
activity.mcmc.median[,i] <- sapply(1:3000, function(idx) 100*(sum(p.alive.mcmc.median[[idx]] >= p.vec[i] & x.star.true[[idx]][,3] >0) + sum(p.alive.mcmc.median[[idx]] < p.vec[i] & x.star.true[[idx]][,3] == 0)) / dataset.characteristics$cohort_size[idx])
activity.mcmc.mean[,i] <- sapply(1:3000, function(idx) 100*(sum(p.alive.mcmc.mean[[idx]] >= p.vec[i] & x.star.true[[idx]][,3] >0) + sum(p.alive.mcmc.mean[[idx]] < p.vec[i] & x.star.true[[idx]][,3] == 0)) / dataset.characteristics$cohort_size[idx])

activity.ind.median[,i] <- sapply(1:3000, function(idx) 100*(sum(p.alive.ind.median[[idx]] >= p.vec[i] & x.star.true[[idx]][,3] >0) + sum(p.alive.ind.median[[idx]] < p.vec[i] & x.star.true[[idx]][,3] == 0)) / dataset.characteristics$cohort_size[idx])
activity.ind.mean[,i] <- sapply(1:3000, function(idx) 100*(sum(p.alive.ind.mean[[idx]] >= p.vec[i] & x.star.true[[idx]][,3] >0) + sum(p.alive.ind.mean[[idx]] < p.vec[i] & x.star.true[[idx]][,3] == 0)) / dataset.characteristics$cohort_size[idx])

activity.MLE[,i] <- sapply(1:3000, function(idx) 100*(sum(p.alive.MLE[[idx]] >= p.vec[i] & x.star.true[[idx]][,3] >0) + sum(p.alive.MLE[[idx]] < p.vec[i] & x.star.true[[idx]][,3] == 0)) / dataset.characteristics$cohort_size[idx])
}

p.max.mcmc.median <- p.vec[which.max(apply(activity.mcmc.median,2,mean))]
p.max.mcmc.mean <- p.vec[which.max(apply(activity.mcmc.mean,2,mean))]
p.max.ind.mean <- p.vec[which.max(apply(activity.ind.mean,2,mean))]
p.max.ind.median <- p.vec[which.max(apply(activity.ind.median,2,mean))]
p.max.MLE <- p.vec[which.max(apply(activity.MLE,2,function(x) mean(x, na.rm=T)))]

par(mfrow=c(1,5))
boxplot(activity.mcmc.median[,which.max(apply(activity.mcmc.median,2,mean))], main= "mcmc median with p.max = 0.55",outline=F,ylim=c(70,100))
boxplot(activity.ind.median[,which.max(apply(activity.ind.median,2,mean))], main= "ind median with p.max = 0.15",outline=F,ylim=c(70,100))
boxplot(activity.mcmc.mean[,which.max(apply(activity.mcmc.mean,2,mean))], main= "mcmc mean with p.max = 0.65",outline=F,ylim=c(70,100))
boxplot(activity.ind.mean[,which.max(apply(activity.ind.mean,2,mean))], main= "ind mean with p.max = 0.3",outline=F,ylim=c(70,100))
boxplot(activity.MLE[,which.max(apply(activity.MLE,2,function(x) mean(x, na.rm=T)))], main= "MLE with p.max = 0.6",outline=F,ylim=c(70,100))


################################################
### step 3: include real data sets
#################################################

load("real_datasets.Rdata")
load("real_est/MLE_est_real.Rdata")
p.alive.mcmc.median.real <- p.alive.ind.median.real <- p.alive.mcmc.mean.real <- p.alive.ind.mean.real <- p.alive.MLE.real <- list()
load("dataset_characteristics_real.Rdata")


P.alive.ind <- function(t.x, T, lambda, mu){
  return((1+(lambda / (lambda + mu)) * (exp((lambda + mu)*(T - t.x)) - 1))^(-1))
}

for (row in 1:22){
  
  purchase_data <- real.datasets[[dataset.characteristics.real$dataset[row]]]
  load(paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("real_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- dataset.characteristics.real$cohort_size[row]
  t <- dataset.characteristics.real$calibration[row]
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  p.alive.mcmc.median.real[[row]] <- apply(MCMC_ind[4,,], 2,median)
  p.alive.mcmc.mean.real[[row]] <- apply(MCMC_ind[4,,], 2,mean)
  
  p.alive.ind.median.real[[row]] <- sapply(1:n, function(x) median(P.alive.ind(cbs[x,2], cbs[x,3], MCMC_ind[1,,x], MCMC_ind[2,,x])))
  p.alive.ind.mean.real[[row]] <- sapply(1:n, function(x) mean(P.alive.ind(cbs[x,2], cbs[x,3], MCMC_ind[1,,x], MCMC_ind[2,,x])))
  
  p.alive.MLE.real[[row]] <- CLVTools:::pnbd_nocov_PAlive(MLE.est[row,1],MLE.est[row,2],MLE.est[row,3],MLE.est[row,4], cbs[,1], cbs[,2], cbs[,3])
  
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop

save(p.alive.ind.mean.real, p.alive.ind.median.real, p.alive.mcmc.mean.real, p.alive.mcmc.median.real, p.alive.MLE.real, file="results/P_alives_real.Rdata")


################################################
### step 4: derive activity status from p(alive) for real data
#################################################

load("results/forecast_true_real.Rdata")
x.star.true.real <- x.star.true

activity.mcmc.median.real <- activity.mcmc.mean.real <- activity.ind.median.real <- activity.ind.mean.real <- activity.MLE.real <- data.frame(rep(NA, 22))

for (i in 1:21){
  activity.mcmc.median.real[,i] <- sapply(1:22, function(idx) 100*(sum(p.alive.mcmc.median.real[[idx]] >= p.vec[i] & x.star.true.real[[idx]][,3] >0) + sum(p.alive.mcmc.median.real[[idx]] < p.vec[i] & x.star.true.real[[idx]][,3] == 0)) / dataset.characteristics.real$cohort_size[idx])
  activity.mcmc.mean.real[,i] <- sapply(1:22, function(idx) 100*(sum(p.alive.mcmc.mean.real[[idx]] >= p.vec[i] & x.star.true.real[[idx]][,3] >0) + sum(p.alive.mcmc.mean.real[[idx]] < p.vec[i] & x.star.true.real[[idx]][,3] == 0)) / dataset.characteristics.real$cohort_size[idx])
  
  activity.ind.median.real[,i] <- sapply(1:22, function(idx) 100*(sum(p.alive.ind.median.real[[idx]] >= p.vec[i] & x.star.true.real[[idx]][,3] >0) + sum(p.alive.ind.median.real[[idx]] < p.vec[i] & x.star.true.real[[idx]][,3] == 0)) / dataset.characteristics.real$cohort_size[idx])
  activity.ind.mean.real[,i] <- sapply(1:22, function(idx) 100*(sum(p.alive.ind.mean.real[[idx]] >= p.vec[i] & x.star.true.real[[idx]][,3] >0) + sum(p.alive.ind.mean.real[[idx]] < p.vec[i] & x.star.true.real[[idx]][,3] == 0)) / dataset.characteristics.real$cohort_size[idx])
  
  activity.MLE.real[,i] <- sapply(1:22, function(idx) 100*(sum(p.alive.MLE.real[[idx]] >= p.vec[i] & x.star.true.real[[idx]][,3] >0) + sum(p.alive.MLE.real[[idx]] < p.vec[i] & x.star.true.real[[idx]][,3] == 0)) / dataset.characteristics.real$cohort_size[idx])
}

p.max.mcmc.median.real <- p.vec[which.max(apply(activity.mcmc.median.real,2,mean))]
p.max.mcmc.mean.real <- p.vec[which.max(apply(activity.mcmc.mean.real,2,mean))]
p.max.ind.mean.real <- p.vec[which.max(apply(activity.ind.mean.real,2,mean))]
p.max.ind.median.real <- p.vec[which.max(apply(activity.ind.median.real,2,mean))]
p.max.MLE.real <- p.vec[which.max(apply(activity.MLE.real,2,function(x) mean(x, na.rm=T)))]


par(mfrow=c(1,5))
boxplot(activity.mcmc.median[,which.max(apply(activity.mcmc.median,2,mean))], main= "mcmc median with p.max = 0.55",outline=F,ylim=c(30,100))
points(rep(1,22),activity.mcmc.median.real[,which.max(apply(activity.mcmc.median.real,2,mean))], pch=19, col="red")

boxplot(activity.ind.median[,which.max(apply(activity.ind.median,2,mean))], main= "ind median with p.max = 0.15",outline=F,ylim=c(30,100))
points(rep(1,22),activity.ind.median.real[,which.max(apply(activity.ind.median.real,2,mean))], pch=19, col="red")

boxplot(activity.mcmc.mean[,which.max(apply(activity.mcmc.mean,2,mean))], main= "mcmc mean with p.max = 0.65",outline=F,ylim=c(30,100))
points(rep(1,22),activity.mcmc.mean.real[,which.max(apply(activity.mcmc.mean.real,2,mean))], pch=19, col="red")

boxplot(activity.ind.mean[,which.max(apply(activity.ind.mean,2,mean))], main= "ind mean with p.max = 0.3",outline=F,ylim=c(30,100))
points(rep(1,22),activity.ind.mean.real[,which.max(apply(activity.ind.mean.real,2,mean))], pch=19, col="red")

boxplot(activity.MLE[,which.max(apply(activity.MLE,2,function(x) mean(x, na.rm=T)))], main= "MLE with p.max = 0.6",outline=F,ylim=c(30,100))
points(rep(1,22),activity.MLE.real[,which.max(apply(activity.MLE.real,2,function(x) mean(x, na.rm=T)))], pch=19, col="red")

