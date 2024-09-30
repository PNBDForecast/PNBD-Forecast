load("results/110_x_star_deviation_real.Rdata")
load("results/forecasts/forecast_true.Rdata")
load("sample_analysis_real.Rdata")
load("results/101_x_star_deviation.Rdata")


par(mfrow=c(1,2))

boxplot(MAPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,2] > 10))),1:7],outline=F, main="MAPE of future purchases (T* = 26)")
for (i in 1:7){
points(rep(i,13),MAPE.real[[2]][sample.analysis.real$dataset !=1,i],col="red",pch=4)
}

boxplot(RMSPE[[2]][unlist(lapply(x.star.true, function(x) sum(x[,2] > 10))),1:7],outline=F, main="RMSPE of future purchases (T* = 26)")
for (i in 1:7){
  points(rep(i,13),RMSPE.real[[2]][sample.analysis.real$dataset !=1,i],col="red",pch=4)
}
