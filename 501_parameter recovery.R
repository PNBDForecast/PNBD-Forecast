library(forecast)
load("results/MLE_est.Rdata")
MLE.est <- as.data.frame(MLE.est)

het.est.median <- het.est.mean <- as.data.frame(matrix(nrow=3000,ncol=4)) 
colnames(het.est.median) <- colnames(het.est.mean) <- c("r","alpha","s","beta")

for (row in 1:3000){
  load(paste("sim_est/MCMC_het/het_est_",row,".Rdata",sep=""))
  het.est.median[row,] <- apply(MCMC_het,2,median)[1:4]
  het.est.mean[row,] <- apply(MCMC_het,2,mean)[1:4]
}


#--------------------------------------------------------------------------------------

# accuracy for r
r.estimates <- cbind(MLE.est$r, het.est.median$r, het.est.mean$r)
colnames(r.estimates) <- c("MLE","median","mean")

r.accuracy <- apply(r.estimates,2,function(x) accuracy(x, sample.analysis$r))
rownames(r.accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(r.accuracy,3)

#--------------------------------------------------------------------------------------

# accuracy for alpha
alpha.estimates <- cbind(MLE.est$alpha, het.est.median$alpha, het.est.mean$alpha)
colnames(alpha.estimates) <- c("MLE","median","mean")

alpha.accuracy <- apply(alpha.estimates,2,function(x) accuracy(x, sample.analysis$alpha))
rownames(alpha.accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(alpha.accuracy,3)

#--------------------------------------------------------------------------------------

# accuracy for s
s.estimates <- cbind(MLE.est$s, het.est.median$s, het.est.mean$s)
colnames(s.estimates) <- c("MLE","median","mean")

s.accuracy <- apply(s.estimates,2,function(x) accuracy(x, sample.analysis$s))
rownames(s.accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(s.accuracy,3)

#--------------------------------------------------------------------------------------

# accuracy for beta
beta.estimates <- cbind(MLE.est$beta, het.est.median$beta, het.est.mean$beta)
colnames(beta.estimates) <- c("MLE","median","mean")

beta.accuracy <- apply(beta.estimates,2,function(x) accuracy(x, sample.analysis$beta))
rownames(beta.accuracy) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(beta.accuracy,3)

#--------------------------------------------------------------------------------------

# accuracy for E(lambda)
E.la.estimates <- cbind(MLE.est$r/MLE.est$alpha, het.est.median$r/het.est.median$alpha, het.est.mean$r/het.est.mean$alpha)
colnames(E.la.estimates) <- c("MLE","median","mean")

E.la.estimates <- apply(E.la.estimates,2,function(x) accuracy(x, sample.analysis$r/sample.analysis$alpha))
rownames(E.la.estimates) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(E.la.estimates,3)


# accuracy for E(mu)
E.mu.estimates <- cbind(MLE.est$s/MLE.est$beta, het.est.median$s/het.est.median$beta, het.est.mean$s/het.est.mean$beta)
colnames(E.mu.estimates) <- c("MLE","median","mean")

E.mu.estimates <- apply(E.mu.estimates,2,function(x) accuracy(x, sample.analysis$s/sample.analysis$beta))
rownames(E.mu.estimates) <- c("ME", "RMSE", "MAE", "MPE", "MAPE")
round(E.mu.estimates,3)

#--------------------------------------------------------------------------------------
round(r.accuracy[c(2,5),],3)
round(alpha.accuracy[c(2,5),],3)
round(s.accuracy[c(2,5),],3)
round(beta.accuracy[c(2,5),],3)


#--------------------------------------------------------------------------------------
# plot s values
plot(c(sample.analysis$s[sample.analysis$s==0.25],sample.analysis$s[sample.analysis$s==1],sample.analysis$s[sample.analysis$s==(32/18)]),type="l",ylim=c(0,3),lwd=2,ylab="",main="s accuracy")
lines(c(MLE.est$s[order(MLE.est$s)][sample.analysis$s==0.25], MLE.est$s[order(MLE.est$s)][sample.analysis$s==1], MLE.est$s[order(MLE.est$s)][sample.analysis$s==(32/18)]),col="blue")
lines(c(het.est.median$s[order(het.est.median$s)][sample.analysis$s==0.25], het.est.median$s[order(het.est.median$s)][sample.analysis$s==1], het.est.median$s[order(het.est.median$s)][sample.analysis$s==(32/18)]),col="red")
abline(v=1000,lwd=4,col="white")
abline(v=2000,lwd=4,col="white")
legend("topleft",col=c(1,"blue","red"),lty=1,legend=c("true", "MLE", "median"))
        
