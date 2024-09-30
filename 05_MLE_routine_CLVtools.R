source("00_general data.R")
load("sample_analysis.Rdata")

library(CLVTools)

MLE.est <- cbind("r" = rep(0,3000), "alpha" = rep(0,3000), 
                 "s" = rep(0,3000), "beta" = rep(0,3000),
                 "LL" = rep(0,3000))



for (row in 1:3000){ 
  load(paste("sim_data/dataset_",row,".Rdata",sep=""))
  
  t <- sample.analysis$calibration[row]
  n <- sample.analysis$cohort.size[row]
  
  transaction.data.frame <- c(0,1)
  
  for (cust in 1:n){
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
  
 save(MLE.est,file="sim_est/MLE_est_CLVtools.Rdata") 
}
