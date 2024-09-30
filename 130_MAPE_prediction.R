load("dataset_characteristics_real.Rdata")
load("results/120_x_star_regression_models.Rdata")
predictors.real <- dataset.characteristics.real[,-1]

real.idx <- dataset.characteristics.real$dataset != 1

# determine MAPE prediction for the real data sets
  
  prediction.sim.13 <- predict(lm.13$lm.sim.13, newdata = predictors.real, interval = "prediction")
  prediction.ind.13 <- predict(lm.13$lm.ind.13, newdata = predictors.real, interval = "prediction")
  prediction.MLE.13 <- predict(lm.13$lm.MLE.13, newdata = predictors.real, interval = "prediction")
  prediction.heur.13 <- predict(lm.13$lm.heur.13, newdata = predictors.real, interval = "prediction")
  
  prediction.sim.26 <- predict(lm.26$lm.sim.26, newdata = predictors.real, interval = "prediction")
  prediction.ind.26 <- predict(lm.26$lm.ind.26, newdata = predictors.real, interval = "prediction")
  prediction.MLE.26 <- predict(lm.26$lm.MLE.26, newdata = predictors.real, interval = "prediction")
  prediction.heur.26 <- predict(lm.26$lm.heur.26, newdata = predictors.real, interval = "prediction")
  
  prediction.sim.52 <- predict(lm.52$lm.sim.52, newdata = predictors.real, interval = "prediction")
  prediction.ind.52 <- predict(lm.52$lm.ind.52, newdata = predictors.real, interval = "prediction")
  prediction.MLE.52 <- predict(lm.52$lm.MLE.52, newdata = predictors.real, interval = "prediction")
  prediction.heur.52 <- predict(lm.52$lm.heur.52, newdata = predictors.real, interval = "prediction")

  
  
# Comparison with true MAPE
# kein signifikanter Unterschiedn zwischen den forecast periods => wähle 
par(mfrow=(c(2,2)))


plot(MAPE.real[[2]]$sim_med[real.idx],type="l",col="red", main="sim (26 weeks)",ylim=c(0,150),xlab="",ylab="")
lines(prediction.sim.26[real.idx,2], lty=2)
lines(prediction.sim.26[real.idx,3], lty=2)
plot(MAPE.real[[2]]$ind_med[real.idx],type="l",col="red", main="ind (26 weeks)",ylim=c(0,150),xlab="",ylab="")
lines(prediction.ind.26[real.idx,2], lty=2)
lines(prediction.ind.26[real.idx,3], lty=2)
plot(MAPE.real[[2]]$MLE[real.idx],type="l",col="red", main="MLE (26 weeks)",ylim=c(0,150),xlab="",ylab="")
lines(prediction.MLE.26[real.idx,2], lty=2)
lines(prediction.MLE.26[real.idx,3], lty=2)
plot(MAPE.real[[2]]$heur[real.idx],type="l",col="red", main="heur (26 weeks)",ylim=c(0,150),xlab="",ylab="")
lines(prediction.heur.26[real.idx,2], lty=2)
lines(prediction.heur.26[real.idx,3], lty=2)
legend("bottom",col=c(1,2),lty=c(2,1),legend=c("prediction interval","true value"), horiz=T)


