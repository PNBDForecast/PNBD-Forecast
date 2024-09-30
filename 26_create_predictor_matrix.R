load("results/x_star_deviations.Rdata")
load("dataset_characteristics.Rdata")
load("Sample_analysis.Rdata")

predictors <- dataset.characteristics[,c(1:2,11:16)]

N <- sum(dataset.characteristics$cohort.size)
load("results/forecast_mcmc.Rdata")
load("results/forecast_true.Rdata")


"cohort.size" <- "calibration" <- "repeat.buyers" <- "mean.F" <- "CV.F" <- "mean.R" <- 
  "CV.R" <- "Gini.x" <- "lambda" <- "mu" <- "z" <- "recency" <- "frequency" <- 0




for (row in 1:3000){
  n <- predictors$cohort.size[row]
  cohort.size <- c(cohort.size, rep(predictors$cohort.size[row],n))
  calibration <- c(calibration, rep(predictors$calibration[row],n))
  repeat.buyers <- c(repeat.buyers, rep(predictors$repeat.buyers[row],n))
  mean.F <- c(mean.F, rep(predictors$mean.F[row],n))
  CV.F <- c(CV.F, rep(predictors$CV.F[row],n))
  mean.R <- c(mean.R, rep(predictors$mean.R[row],n))
  CV.R <- c(CV.R, rep(predictors$CV.R[row],n))
  Gini.x <- c(Gini.x, rep(predictors$Gini.x[row],n))
  
  load(paste("sim_est/MCMC_ind/ind_est_",row,".Rdata",sep="")) # load individual MCMC estimates  
  lambda <- c(lambda, apply(MCMC_ind[1,,],2,median))
  mu <- c(mu, apply(MCMC_ind[2,,],2,median))
  z <- c(z, apply(MCMC_ind[4,,],2,median))
  
  
  load(paste("sim_data/dataset_",row,".Rdata",sep="")) # load purchase data
  n <- sample.analysis$cohort.size[row]
  t <- sample.analysis$calibration[row]
  cbs <- array(NA_real_,dim=c(n,3))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal")
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  recency <- c(recency, cbs[,3] - cbs[,2])
  frequency <- c(frequency, cbs[,1])
  
  cat("end of row",row, format(Sys.time(), '%H:%M:%S'), fill = TRUE)
}

ind.predictors <- data.frame(cohort.size, calibration, repeat.buyers, mean.F, CV.F, mean.R, CV.R, Gini.x, 
                            lambda, mu, z, recency, frequency)[-1,]
rm(purchase_data, calibration, cohort.size, cust, CV.F, CV.R, frequency, Gini.x, lambda, MCMC_ind, mean.F, mean.R, mu, recency, repeat.buyers)

save(predictors, ind.predictors, file="results/predictors.Rdata")


