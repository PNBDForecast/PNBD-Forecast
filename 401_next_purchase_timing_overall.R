library(forecast)
load("2023_sample_analysis.Rdata")

accuracy.median <- accuracy(1:10, 1:10)
accuracy.mean <- accuracy(1:10, 1:10)

t.x1.median.with.mu.total <- t.x1.mean.with.mu.total <- t.x1.median.total <- t.x1.mean.total <- t.x1.true.total <- 0



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
  # latest.tx1[row] <- max(t.x1.true)
  # the median draw for t.x+1 only requires the median draw as it is monotone
  # t.x1.median <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,],2,median))),39)
  t.x1.median.with.mu <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,1:n] + MCMC_ind[2,,1:n],2,median))),39)
  t.x1.median.with.mu <- t.x1.median.with.mu[t.x1.idx]  
  t.x1.mean.with.mu <- 7 * (1/ apply(MCMC_ind[1,,1:n]+MCMC_ind[2,,1:n],2,median))
  t.x1.mean.with.mu <- t.x1.mean.with.mu[t.x1.idx]
  
  t.x1.median <- pmin(7 * ((log(2) / apply(MCMC_ind[1,,1:n] ,2,median))),39)
  t.x1.median <- t.x1.median[t.x1.idx]  
  t.x1.mean <- 7 * (1/ apply(MCMC_ind[1,,1:n],2,median))
  t.x1.mean <- t.x1.mean[t.x1.idx]
  
  t.x1.true.total <- c(t.x1.true.total, t.x1.true)
  t.x1.mean.total <- c(t.x1.mean.total, t.x1.mean)
  t.x1.mean.with.mu.total <- c(t.x1.mean.with.mu.total, t.x1.mean.with.mu)
  t.x1.median.total <- c(t.x1.median.total, t.x1.median)
  t.x1.median.with.mu.total <- c(t.x1.median.with.mu.total, t.x1.median.with.mu)
 
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
} # end of row loop


accuracy.median <- accuracy(t.x1.median[t.x1.true != Inf], t.x1.true[t.x1.true != Inf])
accuracy.mean <- accuracy(t.x1.mean[t.x1.true != Inf], t.x1.true[t.x1.true != Inf])

accuracy.median.with.mu <- accuracy(t.x1.median.with.mu[t.x1.true != Inf], t.x1.true[t.x1.true != Inf])
accuracy.mean.with.mu <- accuracy(t.x1.mean.with.mu[t.x1.true != Inf], t.x1.true[t.x1.true != Inf])

par(mfrow=c(2,4))
boxplot(t.x1.mean.total - t.x1.true.total, main="t.x1.mean signed",ylim=c(-70,70))
boxplot(t.x1.median.total - t.x1.true.total, main="t.x1.median signed",ylim=c(-70,70))
boxplot(t.x1.mean.with.mu.total - t.x1.true.total, main="t.x1.mean with mu signed",ylim=c(-70,70))
boxplot(t.x1.median.with.mu.total - t.x1.true.total, main="t.x1.median with mu signed",ylim=c(-70,70))

boxplot(abs(t.x1.mean.total - t.x1.true.total), main="t.x1.mean unsigned",ylim=c(0,70))
boxplot(abs(t.x1.median.total - t.x1.true.total), main="t.x1.median unsigned",ylim=c(0,70))
boxplot(abs(t.x1.mean.with.mu.total - t.x1.true.total), main="t.x1.mean with mu unsigned",ylim=c(0,70))
boxplot(abs(t.x1.median.with.mu.total - t.x1.true.total), main="t.x1.median with mu unsigned",ylim=c(0,70))

summary(t.x1.mean.total - t.x1.true.total)
summary(t.x1.median.total - t.x1.true.total)
summary(t.x1.mean.with.mu.total - t.x1.true.total)
summary(t.x1.median.with.mu.total - t.x1.true.total)

summary(abs(t.x1.mean.total - t.x1.true.total))
summary(abs(t.x1.median.total - t.x1.true.total))
summary(abs(t.x1.mean.with.mu.total - t.x1.true.total))
summary(abs(t.x1.median.with.mu.total - t.x1.true.total))

boxplot(accuracy.median[,c(1,3)],outline=F, main = "t.x+1 accuracy")
abline(h=0)
