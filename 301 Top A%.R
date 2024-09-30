library(forecast)

load("sample_analysis.Rdata")
load("results/forecasts/forecast_MLE.Rdata")
load("results/forecasts/forecast_MCMC.Rdata")
load("results/forecasts/forecast_heur.Rdata")
load("results/forecasts/forecast_true.Rdata")

# fix T.star to 52 weeks
T.star <- 3

# define helper: vector with active customer ratio
  active.customer.ratio <- rep(NA,3000)
  for (row in 1:3000){
    active.customer.ratio[row] <- 100*(sum(x.star.true[[row]][,3] > 0) / nrow(x.star.true[[row]]))
  }




# Step 1: Identify the true A% best customers

top.10.ncust.true <- top.20.ncust.true <- rep(NA,3000)
top.10.npurch.true <- top.20.npurch.true <- rep(NA,3000)

for (row in 1:3000){
  top.10.npurch.true[row] <- quantile(x.star.true[[row]][,3], 0.9)
  top.20.npurch.true[row] <- quantile(x.star.true[[row]][,3], 0.8)
  
  top.10.ncust.true[row] <- sum(x.star.true[[row]][,3] >= top.10.npurch.true[row]) 
  top.20.ncust.true[row] <- sum(x.star.true[[row]][,3] >= top.20.npurch.true[row]) 
  
}

# Step 2: Identify the predicted A% best customers

  top.10.ncust.sim <- top.10.ncust.ind <- top.10.ncust.MLE <- top.10.ncust.heur <- rep(NA,3000)
  top.20.ncust.sim <- top.20.ncust.ind <- top.20.ncust.MLE <- top.20.ncust.heur <- rep(NA,3000)
  top.10.npurch.sim <- top.10.npurch.ind <- top.10.npurch.MLE <- top.10.npurch.heur <- rep(NA,3000)
  top.20.npurch.sim <- top.20.npurch.ind <- top.20.npurch.MLE <- top.20.npurch.heur <- rep(NA,3000)
  
  for (row in 1:3000){
   
    top.10.npurch.sim[row] <- quantile(x.star.mcmc[[row]][3,,1,2], 0.9)
    top.10.npurch.ind[row] <- quantile(x.star.mcmc[[row]][3,,2,2], 0.9)
    top.10.npurch.MLE[row] <- quantile(x.star.MLE[[row]][,3], 0.9, na.rm=T)
    top.10.npurch.heur[row] <- quantile(x.star.heur[[row]][,3], 0.9)
    
    top.20.npurch.sim[row] <- quantile(x.star.mcmc[[row]][3,,1,2], 0.8)
    top.20.npurch.ind[row] <- quantile(x.star.mcmc[[row]][3,,2,2], 0.8)
    top.20.npurch.MLE[row] <- quantile(x.star.MLE[[row]][,3], 0.8, na.rm=T)
    top.20.npurch.heur[row] <- quantile(x.star.heur[[row]][,3], 0.8)
    
    top.10.ncust.sim[row] <- sum(x.star.mcmc[[row]][3,,1,2] >= top.10.npurch.sim[row])
    top.10.ncust.ind[row] <- sum(x.star.mcmc[[row]][3,,2,2] >= top.10.npurch.ind[row])
    top.10.ncust.MLE[row] <- sum(x.star.MLE[[row]][,3] >= top.10.npurch.MLE[row], na.rm=T)
    top.10.ncust.heur[row] <- sum(x.star.heur[[row]][,3] >= top.10.npurch.heur[row])
    
    top.20.ncust.sim[row] <- sum(x.star.mcmc[[row]][3,,1,2] >= top.20.npurch.sim[row]) 
    top.20.ncust.ind[row] <- sum(x.star.mcmc[[row]][3,,2,2] >= top.20.npurch.ind[row]) 
    top.20.ncust.MLE[row] <- sum(x.star.MLE[[row]][,3] >= top.20.npurch.MLE[row], na.rm=T) 
    top.20.ncust.heur[row] <- sum(x.star.heur[[row]][,3] >= top.20.npurch.heur[row]) 
    
  }

# Step 3: Compare true and predicted A% customers
  
  correct.top.10.sim <- correct.top.10.ind <- correct.top.10.MLE <- correct.top.10.heur <- rep(NA, 3000)
  correct.top.20.sim <- correct.top.20.ind <- correct.top.20.MLE <- correct.top.20.heur <- rep(NA, 3000)
  
  for (row in 1:3000){
   
    correct.top.10.ind[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true[row]], 
                                                order(x.star.mcmc[[row]][3,,2,2], decreasing = T)[1:top.10.ncust.ind[row]])) / top.10.ncust.true[row]
    
    if (!is.na(x.star.MLE[[row]][1,3])) correct.top.10.MLE[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true[row]], 
                                                    order(x.star.MLE[[row]][,3], decreasing = T)[1:top.10.ncust.MLE[row]])) / top.10.ncust.true[row]
    
    correct.top.10.heur[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.10.ncust.true[row]], 
                                                    order(x.star.heur[[row]][,3], decreasing = T)[1:top.10.ncust.heur[row]])) / top.10.ncust.true[row]
    ###########################
    
    correct.top.20.ind[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true[row]], 
                                                    order(x.star.mcmc[[row]][3,,2,2], decreasing = T)[1:top.20.ncust.ind[row]])) / top.20.ncust.true[row]
    
    if (!is.na(x.star.MLE[[row]][1,3])) correct.top.20.MLE[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true[row]], 
                                                    order(x.star.MLE[[row]][,3], decreasing = T)[1:top.20.ncust.MLE[row]])) / top.20.ncust.true[row]
    
    correct.top.20.heur[row] <- 100*length(intersect(order(x.star.true[[row]][,3], decreasing = T)[1:top.20.ncust.true[row]], 
                                                     order(x.star.heur[[row]][,3], decreasing = T)[1:top.20.ncust.heur[row]])) / top.20.ncust.true[row]
    
    
  }
  
rm(x.star.heur, x.star.mcmc, x.star.mcmc, x.star.true, row, T.star)
save.image("results/301_Top_A%.Rdata")  
