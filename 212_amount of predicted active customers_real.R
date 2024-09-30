amt.true <- unlist(lapply(x.star.true, function(x) sum(x[,2] > 0)))
amt.sim <- unlist(lapply(x.star.MCMC, function(x) sum(x[2,,1,2] >= 0.5)))
amt.ind <- unlist(lapply(x.star.MCMC, function(x) sum(x[2,,2,2] >= 0.5)))
amt.MLE <- unlist(lapply(x.star.MLE, function(x) sum(x[,2] >= 0.5)))
amt.heur <- unlist(lapply(x.star.heur, function(x) sum(x[,2] >= 0.5)))
                   
amt.active.cust <- data.frame(sample.analysis.real, amt.true, amt.sim, amt.ind, amt.MLE, amt.heur)
