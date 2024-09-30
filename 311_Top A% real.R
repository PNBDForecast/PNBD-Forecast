
load("sample_analysis_real.Rdata")
load("results/forecasts/forecast_MLE_real.Rdata")
load("results/forecasts/forecast_MCMC_real.Rdata")
load("results/forecasts/forecast_heur_real.Rdata")
load("results/forecasts/forecast_true_real.Rdata")

# fix T.star to 52 weeks
T.star <- 3



# Step 1: Identify the true A% best customers

top.10.ncust.true.real <- top.20.ncust.true.real <- rep(NA,22)
top.10.npurch.true.real <- top.20.npurch.true.real <- rep(NA,22)

for (row in 1:22){
  top.10.npurch.true.real[row] <- quantile(x.star.true[[row]][,3], 0.9)
  top.20.npurch.true.real[row] <- quantile(x.star.true[[row]][,3], 0.8)
  
  top.10.ncust.true.real[row] <- sum(x.star.true[[row]][,3] >= top.10.npurch.true.real[row]) 
  top.20.ncust.true.real[row] <- sum(x.star.true[[row]][,3] >= top.20.npurch.true.real[row]) 
  
}

# Step 2: Identify the predicted A% best customers

top.10.ncust.sim.real <- top.10.ncust.ind.real <- top.10.ncust.MLE.real <- top.10.ncust.heur.real <- rep(NA,22)
top.20.ncust.sim.real <- top.20.ncust.ind.real <- top.20.ncust.MLE.real <- top.20.ncust.heur.real <- rep(NA,22)
top.10.npurch.sim.real <- top.10.npurch.ind.real <- top.10.npurch.MLE.real <- top.10.npurch.heur.real <- rep(NA,22)
top.20.npurch.sim.real <- top.20.npurch.ind.real <- top.20.npurch.MLE.real <- top.20.npurch.heur.real <- rep(NA,22)

for (row in 1:22){
  
  top.10.npurch.sim.real[row] <- quantile(x.star.MCMC[[row]][3,,1,2], 0.9)
  top.10.npurch.ind.real[row] <- quantile(x.star.MCMC[[row]][3,,2,2], 0.9)
  top.10.npurch.MLE.real[row] <- quantile(x.star.MLE[[row]][,3], 0.9, na.rm=T)
  top.10.npurch.heur.real[row] <- quantile(x.star.heur[[row]][,3], 0.9)
  
  top.20.npurch.sim.real[row] <- quantile(x.star.MCMC[[row]][3,,1,2], 0.8)
  top.20.npurch.ind.real[row] <- quantile(x.star.MCMC[[row]][3,,2,2], 0.8)
  top.20.npurch.MLE.real[row] <- quantile(x.star.MLE[[row]][,3], 0.8, na.rm=T)
  top.20.npurch.heur.real[row] <- quantile(x.star.heur[[row]][,3], 0.8)
  
  top.10.ncust.sim.real[row] <- sum(x.star.MCMC[[row]][3,,1,2] >= top.10.npurch.sim.real[row])
  top.10.ncust.ind.real[row] <- sum(x.star.MCMC[[row]][3,,2,2] >= top.10.npurch.ind.real[row])
  top.10.ncust.MLE.real[row] <- sum(x.star.MLE[[row]][,3] >= top.10.npurch.MLE.real[row], na.rm=T)
  top.10.ncust.heur.real[row] <- sum(x.star.heur[[row]][,3] >= top.10.npurch.heur.real[row])
  
  top.20.ncust.sim.real[row] <- sum(x.star.MCMC[[row]][3,,1,2] >= top.20.npurch.sim.real[row]) 
  top.20.ncust.ind.real[row] <- sum(x.star.MCMC[[row]][3,,2,2] >= top.20.npurch.ind.real[row]) 
  top.20.ncust.MLE.real[row] <- sum(x.star.MLE[[row]][,3] >= top.20.npurch.MLE.real[row], na.rm=T) 
  top.20.ncust.heur.real[row] <- sum(x.star.heur[[row]][,3] >= top.20.npurch.heur.real[row]) 
  
}

# Step 3: Compare true and predicted A% customers

correct.top.10.sim.real <- correct.top.10.ind.real <- correct.top.10.MLE.real <- correct.top.10.heur.real <- rep(NA, 22)
correct.top.20.sim.real <- correct.top.20.ind.real <- correct.top.20.MLE.real <- correct.top.20.heur.real <- rep(NA, 22)

for (row in 1:22){
  
  correct.top.10.ind.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true.real[row]], 
                                                  order(x.star.MCMC[[row]][3,,2,2], decreasing = T)[1:top.10.ncust.ind.real[row]])) / top.10.ncust.true.real[row]
  
  if (!is.na(x.star.MLE[[row]][1,3])) correct.top.10.MLE.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true.real[row]], 
                                                                                      order(x.star.MLE[[row]][,3], decreasing = T)[1:top.10.ncust.MLE.real[row]])) / top.10.ncust.true.real[row]
  
  correct.top.10.heur.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true.real[row]], 
                                                   order(x.star.heur[[row]][,3], decreasing = T)[1:top.10.ncust.heur.real[row]])) / top.10.ncust.true.real[row]
  ###########################
  
  correct.top.20.ind.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true.real[row]], 
                                                  order(x.star.MCMC[[row]][3,,2,2], decreasing = T)[1:top.20.ncust.ind.real[row]])) / top.20.ncust.true.real[row]
  
  if (!is.na(x.star.MLE[[row]][1,3])) correct.top.20.MLE.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true.real[row]], 
                                                                                      order(x.star.MLE[[row]][,3], decreasing = T)[1:top.20.ncust.MLE.real[row]])) / top.20.ncust.true.real[row]
  
  correct.top.20.heur.real[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true.real[row]], 
                                                   order(x.star.heur[[row]][,3], decreasing = T)[1:top.20.ncust.heur.real[row]])) / top.20.ncust.true.real[row]
  
  
}

rm(x.star.heur, x.star.MCMC, x.star.MLE, x.star.true, row, T.star)
save.image("results/311_Top_A%_real.Rdata")  

