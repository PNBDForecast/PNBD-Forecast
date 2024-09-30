load("results/211_active_customers_optimised_p_real.Rdata")
load("results/212_active_customers_real.Rdata")
load("sample_analysis_real.Rdata")


par(mfrow=c(2,2))
plot(tot.correct.real[[2]]$sim_med[sample.analysis.real$dataset !=1], main="sim for T* = 26 weeks",type="l",ylab="",xlab="",ylim=c(70,90))
lines((opt.correct.real[[2]]$sim.med[sample.analysis.real$dataset !=1]), col="red")

plot(tot.correct.real[[2]]$ind_med[sample.analysis.real$dataset !=1], main="ind for T* = 26 weeks",type="l",ylab="",xlab="",ylim=c(70,90))
lines(opt.correct.real[[2]]$ind.med[sample.analysis.real$dataset !=1], col="red")

plot(tot.correct.real[[2]]$MLE[sample.analysis.real$dataset !=1], main="MLE for T* = 26 weeks",type="l",ylab="",xlab="",ylim=c(70,90))
lines(opt.correct.real[[2]]$MLE[sample.analysis.real$dataset !=1], col="red")

plot(tot.correct.real[[2]]$heur[sample.analysis.real$dataset !=1], main="heur for T* = 26 weeks", type="l", ylab="",xlab="",ylim=c(70,90))
lines(opt.correct.real[[2]]$heur[sample.analysis.real$dataset !=1], col="red")