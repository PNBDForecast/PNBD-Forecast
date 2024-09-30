load("results/forecasts/forecast_true.Rdata")
load("sample_analysis.Rdata")
load("results/22_forecast_mcmc_quantiles.Rdata")


fc.period <- c(13,26,52)


correct.interquartile.ratio <- data.frame("75%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  correct.interquartile.ratio[row,1] <- 100*(sum(x.star.true[[row]][,2] >= fc.mcmc[[row]][2,,2] & x.star.true[[row]][,2] <= fc.mcmc[[row]][2,,4]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio[row,2] <- 100*(sum(x.star.true[[row]][,2] >= fc.mcmc[[row]][2,,1] & x.star.true[[row]][,2] <= fc.mcmc[[row]][2,,5]) ) / nrow(x.star.true[[row]])
  correct.interquartile.ratio[row,3] <- 100*(sum(x.star.true[[row]][,2] >= (fc.mcmc[[row]][2,,6]-2*fc.mcmc[[row]][2,,7]) & x.star.true[[row]][,2] <= (fc.mcmc[[row]][2,,6]+2*fc.mcmc[[row]][2,,7])) ) / nrow(x.star.true[[row]])
  
}

correct.interquartile.ratio.active.only <- data.frame("75%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  correct.interquartile.ratio.active.only[row,1] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,2] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,4]) ) / sum(x.star.true[[row]][,2] > 0)
  correct.interquartile.ratio.active.only[row,2] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,1] & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,5]) ) / sum(x.star.true[[row]][,2] > 0)
  correct.interquartile.ratio.active.only[row,3] <- 100*(sum(x.star.true[[row]][x.star.true[[row]][,2] > 0,2] >= (fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,6]-2*fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7]) & x.star.true[[row]][x.star.true[[row]][,2] > 0,2] <= (fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,6]+2*fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7])) ) / sum(x.star.true[[row]][,2] > 0)
  
}

# determine relative range to judge on the interval width
relative.range.active <- data.frame("75%" = rep(NA,3000),"90%" = rep(NA,3000), "mean+2sd" = rep(NA,3000))
for (row in 1:3000){
  # x.star.true[[row]][,2] > 0
  relative.range.active[row,1] <- mean((fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,4] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.range.active[row,2] <- mean((fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,5] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,1]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.range.active[row,3] <- mean((4* fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  
}


relative.active.range.75 <- 0
relative.active.range.90 <- 0
relative.active.range.sd <- 0

for (row in 1:3000){
  # x.star.true[[row]][,2] > 0
  relative.active.range.75 <- c(relative.active.range.75,(fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,4] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,2]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.active.range.90 <- c(relative.active.range.90,(fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,5] - fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,1]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
  relative.active.range.sd <- c(relative.active.range.sd,(4* fc.mcmc[[row]][2,x.star.true[[row]][,2] > 0,7]) / x.star.true[[row]][x.star.true[[row]][,2] > 0,2])
}
relative.active.range.75 <- relative.active.range.75[-1]
relative.active.range.90 <- relative.active.range.90[-1]
relative.active.range.sd <- relative.active.range.sd[-1]

save.image("results/103_range_accuracy_mcmc.Rdata")
