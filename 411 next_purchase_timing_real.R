library(forecast)
load("sample_analysis_real.Rdata")
load("real_datasets.Rdata")

accuracy.median <- accuracy(1:10, 1:10)


latest.tx1 <- vector(length=22)

for (row in 1:22){
  
  purchase_data <- real.datasets[[sample.analysis.real$dataset[row]]]
  load(paste("real_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates     
  load(paste("real_est/MCMC_het/het_est_",row,".Rdata",sep="")) # load heterogeneous MCMC estimates
  
  n <- sample.analysis.real$cohort_size[row]
  t <- sample.analysis.real$calibration[row]
  

  # create cbs
  cbs <- array(NA_real_,dim=c(n,3))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal")
  
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  
  t.x1.true <- unlist(lapply(purchase_data, function(x) min(x[x > (7*t) & x < 7*(t+130)])) )
  t.x1.idx <- t.x1.true !=Inf
  t.x1.true <- t.x1.true[t.x1.idx] - (7*t)
  latest.tx1[row] <- max(t.x1.true)
  # the median draw for t.x+1 only requires the median draw as it is monotone
  t.x1.median <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,],2,median))),130)
  t.x1.median <- t.x1.median[t.x1.idx]  
  t.x1.mean <- 7 * (1/ apply(MCMC_ind[1,,],2,median))
  t.x1.mean <- t.x1.mean[t.x1.idx]
  
  accuracy.median <- rbind(accuracy(t.x1.median[t.x1.true != Inf], t.x1.true[t.x1.true != Inf]), accuracy.median)
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop



accuracy.median <- as.data.frame(accuracy.median[1:22,])
# accuracy.mean <- as.data.frame(accuracy.mean[1:22,])


save (accuracy.median, latest.tx1, file="results/accuracy_t_x1_real.Rdata")

accuracy.median.real <- accuracy.median
latest.tx1.real <- latest.tx1
load("results/accuracy_t_x1.Rdata")



#### Plots
par(mfrow=c(1,2))
boxplot(accuracy.median$MAE, main="MAE of t.x+1 accuracy",col="lightpink",ylim=c(0,180), outline=F)
points(rep(1,13),accuracy.median.real$MAE[sample.analysis.real$dataset !=1],col="red",pch=19)
boxplot(latest.tx1, main="latest t.x+1 per data set")
points(rep(1,13), latest.tx1.real[sample.analysis.real$dataset !=1], pch=19, col="blue")

