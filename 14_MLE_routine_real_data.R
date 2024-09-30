source("00_general data.R")
load("sample_analysis_real.Rdata")
load("real_datasets.Rdata")

library(CLVTools)

MLE.est <- cbind("r" = rep(0,22), "alpha" = rep(0,22), 
                 "s" = rep(0,22), "beta" = rep(0,22),
                 "LL" = rep(0,22))



for (row in 1:22){ 
  purchase_data <- real.datasets[[sample.analysis.real$dataset[row]]]
  
  t <- sample.analysis.real$calibration[row]
  n <- sample.analysis.real$cohort_size[row]
  
  purchase_data <- purchase_data[1:n]
  
  transaction.data.frame <- c(0,1)
  
  for (cust in 1:length(purchase_data)){
    TA.temp <- as.data.frame(cbind(rep(cust,times=length(purchase_data[[cust]])),purchase_data[[cust]]))
    transaction.data.frame <- rbind(transaction.data.frame,TA.temp)
  }
  
  transaction.data.frame <- as.data.frame(transaction.data.frame[-1,])
  colnames(transaction.data.frame) <- c("Id","Date")
  transaction.data.frame <- cbind(transaction.data.frame,"Price"=0)
  transaction.data.frame$Date <- round(transaction.data.frame$Date)
  transaction.data.frame$Date <- as.Date(transaction.data.frame$Date,origin = "2015-01-01")
  
  clv.data <- clvdata(transaction.data.frame,time.unit="weeks",date.format="ymd",estimation.split = t)
  
  MLE.est.temp <- pnbd(clv.data)
  MLE.est[row,1:4] <- summary(MLE.est.temp)$coefficients[,1]
  MLE.est[row,5] <- summary(MLE.est.temp)$estimated.LL
  
  if(is.na(MLE.est)[row,5]) {
    cat("try again with Nelder-Mead", fill = TRUE)
    clv.data <- clvdata(transaction.data.frame,time.unit="weeks",date.format="ymd",estimation.split = t)
    
    MLE.est.temp <- pnbd(clv.data,optimx.args=list(method="Nelder-Mead"))
    MLE.est[row,1:4] <- summary(MLE.est.temp)$coefficients[,1]
    MLE.est[row,5] <- summary(MLE.est.temp)$estimated.LL
    
  }
  
  cat("row",row,format(Sys.time(), '%H:%M:%S'), fill = TRUE)
  
 save(MLE.est,file="real_est/MLE_est_real.Rdata") 
}
