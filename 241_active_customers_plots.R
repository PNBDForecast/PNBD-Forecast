# Step 1: use data sets with a repurchaser ratio > 5%

load("results/forecast_heur.Rdata")
load("results/forecast_MCMC.Rdata")
load("results/forecast_MLE.Rdata")
load("results/forecast_true.Rdata")
load("results/rebuyer_sensitivity_precision.Rdata")

rebuyer.ratio <- vector()
for (row in 1:3000){
  rebuyer.ratio[row] <-   100*sum(x.star.true[[row]][,3] > 0)/nrow(x.star.true[[row]])
}

sum(rebuyer.ratio < 1)
# 343



sensitivity.adj <- sensitivity[[3]][which(rebuyer.ratio >= 1),]
precision.adj <- precision[[3]][which(rebuyer.ratio >= 1),]


plot(sensitivity.adj$heur, precision.adj$heur, col="gray31",ylim=c(0,1), xlim=c(0,1), main="precision vs. sensitivity of repurchaser identification",xlab="",ylab="")
points(sensitivity.adj$MLE, precision.adj$MLE, col="blue")
# points(sensitivity.adj$ind_mean, precision.adj$ind_mean, col="green")
points(sensitivity.adj$ind_med, precision.adj$ind_med, col="red")
# points(sensitivity.adj$sim_med, precision.adj$sim_med, col="purple")
legend("bottomleft",legend=c("ind.median","MLE","heur"),col=c("red","blue","gray31"),pch = 1,bty="n")


plot(sort(sensitivity.adj$sim_med + precision.adj$sim_med),type="l",col="purple", main="sensitivity + precision of repurchaser identification",lwd=1.5,ylab="", xlab="")
# lines(sort(sensitivity.adj$sim_mean + precision.adj$sim_mean),type="l",col="red",lty=2,lwd=1.5)
# lines(sort(sensitivity.adj$ind_med + precision.adj$ind_med),type="l",col="purple",lty=1)
# lines(sort(sensitivity.adj$ind_mean + precision.adj$ind_mean),type="l",col="purple",lty=2)
# lines(sort(sensitivity.adj$het_med + precision.adj$het_med),type="l",col="green",lty=2)
# lines(sort(sensitivity.adj$het_mean + precision.adj$het_mean),type="l",col="green",lty=2)
lines(sort(sensitivity.adj$ind_med + precision.adj$ind_med),type="l",col="red",lty=1)
lines(sort(sensitivity.adj$MLE + precision.adj$MLE),type="l",col="blue",lty=1,lwd=1.5)
legend("bottomright",legend=c("ind.med","sim.med","MLE","heur"),col=c("red","purple","blue","gray31"), lty = 1)
lines(sort(sensitivity.adj$heur + precision.adj$heur),type="l",col="gray31",lty=1,lwd=1.5)

par(mfrow=c(2,3))
boxplot(sensitivity.adj$ind_med, main="repeat buyer sensitivity for ind.median",ylim=c(0,1) )
boxplot(sensitivity.adj$MLE, main="repeat buyer sensitivity for MLE",ylim=c(0,1) )
boxplot(sensitivity.adj$heur, main="repeat buyer sensitivity for heuristic",ylim=c(0,1) )
boxplot(precision.adj$ind_med, main="repeat buyer precision for ind.mean",ylim=c(0,1) )
boxplot(precision.adj$MLE, main="repeat buyer precision for MLE",ylim=c(0,1) )
boxplot(precision.adj$heur, main="repeat buyer precision for heuristic",ylim=c(0,1) )