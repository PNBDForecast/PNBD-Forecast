library(BTYD)

# draw_future_transactions <- function(cal.cbs,draws,T.star=cal.cbs$T.star,sample_size=NULL)
# pnbd.ConditionalExpectedTransactions(params,T.star,x,t.x,T.cal,hardie = TRUE)


##############################
# FIXED PARAMETERS
##############################

fc.period <- c(13,26,52)

##############################
# HOUSEKEEPING
##############################

load("sample_analysis.Rdata")
load("sim_est/MLE_est_CLVTools.Rdata")
x.star.MLE <- list()


for (row in 1:3000){
  
  ##########################################
  # PREPARE DATA SET
  ##########################################
  
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
  
  
    
  ##########################################
  # het closed form forecast
  ##########################################
  
  x.star.MLE[[row]] <- matrix(nrow=n, ncol=3)
  colnames(x.star.MLE[[row]]) <- fc.period
  
  for (T.star in 1:3){
    if (!is.na(MLE.est[row,1])){
            x.star.MLE[[row]][,T.star] <- pnbd.ConditionalExpectedTransactions(MLE.est[row,1:4],T.star=fc.period[T.star],cbs[,1],cbs[,2],cbs[,3],hardie = TRUE)
      }
    }

  cat(format(Sys.time(), '%H:%M:%S'),"Datensatz",row,"beendet.","\n")
  
} # end of row

save(x.star.MLE,file="results/forecast_MLE.Rdata")
