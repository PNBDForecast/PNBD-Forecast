
load("sample_analysis_real.Rdata")
load("results/forecasts/forecast_MLE_real.Rdata")
load("results/forecasts/forecast_MCMC_real.Rdata")
load("results/forecasts/forecast_heur_real.Rdata")
load("results/forecasts/forecast_true_real.Rdata")


amt.true <- unlist(lapply(x.star.true, function(x) sum(x[,2] > 0)))
amt.sim <- unlist(lapply(x.star.MCMC, function(x) sum(x[2,,1,2] >= 0.5)))
amt.ind <- unlist(lapply(x.star.MCMC, function(x) sum(x[2,,2,2] >= 0.5)))
amt.MLE <- unlist(lapply(x.star.MLE, function(x) sum(x[,2] >= 0.5)))
amt.heur <- unlist(lapply(x.star.heur, function(x) sum(x[,2] >= 0.5)))

amt.active.cust <- data.frame(sample.analysis.real, amt.true, amt.sim, amt.ind, amt.MLE, amt.heur)

plot((100*amt.sim/amt.true-100)[sample.analysis.real$dataset !=1], type="l", col="red", main="sum of predicted active customers",ylab="",xlab="",ylim=c(-100,100))
lines((100*amt.ind/amt.true-100)[sample.analysis.real$dataset !=1],col="forestgreen")
lines((100*amt.MLE/amt.true-100)[sample.analysis.real$dataset !=1], col="blue")
lines((100*amt.heur/amt.true-100)[sample.analysis.real$dataset !=1], col="gray43")
abline(h=0,lty=2)