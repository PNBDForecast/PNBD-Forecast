library(forecast)
load("2023_sample_analysis.Rdata")

accuracy.median <- accuracy(1:10, 1:10)
accuracy.mean <- accuracy(1:10, 1:10)

latest.tx1 <- vector(length=3000)

for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- sample.analysis$cohort_size[row]
  t <- sample.analysis$calibration[row]
  

  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  
  t.x1.true <- unlist(lapply(purchase_data, function(x) min(x[x > (7*t) & x < 7*(t+39)])) )
  t.x1.idx <- t.x1.true !=Inf
  t.x1.true <- t.x1.true[t.x1.idx] - (7*t)
  latest.tx1[row] <- max(t.x1.true)
  # the median draw for t.x+1 only requires the median draw as it is monotone
  # t.x1.median <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,],2,median))),39)
  t.x1.median.with.mu <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,] + MCMC_ind[2,,],2,median))),39)
  t.x1.median.with.mu <- t.x1.median.with.mu[t.x1.idx]  
  t.x1.mean.with.mu <- 7 * (1/ apply(MCMC_ind[1,,]+MCMC_ind[2,,],2,median))
  t.x1.mean.with.mu <- t.x1.mean.with.mu[t.x1.idx]
  
  t.x1.median <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,] ,2,median))),39)
  t.x1.median <- t.x1.median[t.x1.idx]  
  t.x1.mean <- 7 * (1/ apply(MCMC_ind[1,,],2,median))
  t.x1.mean <- t.x1.mean[t.x1.idx]
  
  accuracy.median <- rbind(accuracy(t.x1.median[t.x1.true != Inf], t.x1.true[t.x1.true != Inf]), accuracy.median)
  accuracy.mean <- rbind(accuracy(t.x1.mean[t.x1.true != Inf], t.x1.true[t.x1.true != Inf]), accuracy.mean)
  
  accuracy.median.with.mu <- rbind(accuracy(t.x1.median.with.mu[t.x1.true != Inf], t.x1.true[t.x1.true != Inf]), accuracy.median)
  accuracy.mean.with.mu <- rbind(accuracy(t.x1.mean.with.mu[t.x1.true != Inf], t.x1.true[t.x1.true != Inf]), accuracy.mean)
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop



accuracy.median <- as.data.frame(accuracy.median[1:3000,])
accuracy.mean <- as.data.frame(accuracy.mean[1:3000,])


accuracy.median.with.mu <- as.data.frame(accuracy.median.with.mu[1:3000,])
accuracy.mean.with.mu <- as.data.frame(accuracy.mean.with.mu[1:3000,])


save (accuracy.mean, accuracy.median,accuracy.mean.with.mu,accuracy.median.with.mu, file="results/accuracy_t_x1.Rdata")

summary(accuracy.median$RMSE)
summary(accuracy.median$MAPE)

summary(accuracy.mean$RMSE)
summary(accuracy.mean$MAPE)


#### Plots
boxplot(accuracy.median$MAE, main="MAE of t.x+1 accuracy")
boxplot(latest.tx1, main="latest t.x+1 per data set")

lm <- lm(accuracy.median$MAE ~ sample.analysis$calibration + sample.analysis$cohort_size + sample.analysis$PR + sample.analysis$DR + sample.analysis$CV_lambda + sample.analysis$CV_mu)
summary(lm)
