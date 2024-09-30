load("sample_analysis.Rdata")
source ("00_general data.R")


for (row in 1:3000){
  r <- sample.analysis$r[row]
  a <- sample.analysis$alpha[row]
  s <- sample.analysis$s[row]
  b <- sample.analysis$beta[row]
  n <- sample.analysis$cohort.size[row]
  
  cat("row",row,"\n")
  
  individual_data <- array(NA_real_,dim=c(n,3))
  dimnames(individual_data)[[2]] <- c("lambda","mu","tau")
  
  #draw lambda
  individual_data[,1] <- rgamma(n,shape=r,rate=a)
  
  #draw mu and tau
  individual_data[,2] <- rgamma(n,shape=s,rate=b)
  individual_data[,3] <- rexp(n,individual_data[,2])
  
  #create purchase data with unit "day", initial purchase in [0,6]
  purchase_data <- list()
  for(cust in 1:n){
    lambda <- individual_data[cust,1]
    mu     <- individual_data[cust,2]
    tau    <- individual_data[cust,3]
    cust_pur <- round(runif(1,-0.5,90.5),0) # generate initial purchase
    while(tail(cust_pur,1)<= min(cust_pur[1]+7*tau,365*4)){ # purchase data for four years
      cust_pur <- c(cust_pur,tail(cust_pur,1)+round(7*rexp(1,lambda),0))
    } #end of while loop
    purchase_data[[cust]] <- cust_pur[1:(length(cust_pur)-1)]
  } # end of cust loop 
  
  cat(format(Sys.time(), '%H:%M:%S'),"row", row, "finished",'\n')
  
  save(individual_data,purchase_data,file=paste("sim_data/dataset_",row,".Rdata",sep=""))
}


          
          
          
          
          
        
   
