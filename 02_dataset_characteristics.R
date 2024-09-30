#####################################################
## DATASET CHARACTERISTICS FOR SIMULATED DATA SETS
#####################################################
library(DescTools)

load("sample_analysis.Rdata")

Gini.x <- repeat.buyers <- mean.F <- CV.F <- mean.R <- CV.R <- rep(0,3000)

for (row in 1:3000){
  
  load(paste("sim_data/dataset_",row,".Rdata",sep=""))
  n <- sample.analysis$cohort.size[row]
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
  
  repeat.buyers[row] <- 100*(sum(cbs[,1] > 0) / n)
  mean.F[row] <- mean(cbs[,1]/cbs[,3])
  CV.F[row] <- sd(cbs[,1])/mean(cbs[,1])
  mean.R[row] <- mean(cbs[,3] - cbs[,2])
  CV.R[row] <- sd(cbs[,3] - cbs[,2]) / mean(cbs[,3] - cbs[,2])
  Gini.x[row] <- Gini(cbs[,1])
}

dataset.characteristics <- data.frame(sample.analysis, repeat.buyers, mean.F, CV.F, mean.R, CV.R, Gini.x)

save(dataset.characteristics, file="dataset_characteristics.Rdata")


#####################################################
## DATASET CHARACTERISTICS FOR REAL DATA SETS
#####################################################

load("sample_analysis_real.Rdata")
load("real_datasets.Rdata")


Gini.x <- repeat.buyers <- mean.F <- CV.F <- mean.R <- CV.R <- rep(0,22)

for (row in 1:22){
  
  purchase_data <- real.datasets[[sample.analysis.real$dataset[row]]]
  
  n <- sample.analysis.real$cohort.size[row]
  t <- sample.analysis.real$calibration[row]
  
  # create cbs
  cbs <- array(NA_real_,dim=c(n,6))
  dimnames(cbs)[[2]] <- c("x","t.x","T.cal","lambda","mu","tau")
  
  cbs[,4:6] <- individual_data[1:n,1:3]
  for (cust in 1:n){
    cbs[cust,1] <- length(unique(purchase_data[[cust]][purchase_data[[cust]]<=7*t]))-1
    cbs[cust,2] <- (max(purchase_data[[cust]][purchase_data[[cust]]<=7*t])-purchase_data[[cust]][1])/7
    cbs[cust,3] <- t - purchase_data[[cust]][1]/7
  } # end of cust loop
  
  repeat.buyers[row] <- 100*(sum(cbs[,1] > 0) / n)
  mean.F[row] <- mean(cbs[,1] / cbs[,3])
  CV.F[row] <- sd(cbs[,1])/mean(cbs[,1])
  mean.R[row] <- mean(cbs[,3] - cbs[,2])
  CV.R[row] <- sd(cbs[,3] - cbs[,2]) / mean(cbs[,3] - cbs[,2])
  Gini.x[row] <- Gini(cbs[,1])
}

dataset.characteristics.real <- data.frame(sample.analysis.real, repeat.buyers, mean.F, CV.F, mean.R, CV.R, Gini.x)

save(dataset.characteristics.real, file="dataset_characteristics_real.Rdata")

