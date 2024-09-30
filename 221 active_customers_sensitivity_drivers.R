load("results/201_active_customers.Rdata")
load("dataset_characteristics.Rdata")


predictors.het <- dataset.characteristics[,c(1:2,11:16)]

#---- 13 weeks-----------

lm.sim.13 <- lm(tot.correct[[1]]$sim_med ~ ., data= predictors.het)
lm.ind.13 <- lm(tot.correct[[1]]$ind_med ~ ., data = predictors.het)
lm.MLE.13 <- lm(tot.correct[[1]]$MLE ~ ., data = predictors.het)
lm.heur.13 <- lm(tot.correct[[1]]$heur ~., data=predictors.het)

drivers.13 <- round(data.frame("est.sim" = lm.sim.13$coefficients[-1], "p.sim" = summary(lm.sim.13)[[4]][-1,4],
                               "est.ind" = lm.ind.13$coefficients[-1], "p.ind" = summary(lm.ind.13)[[4]][-1,4],
                               "est.MLE" = lm.ind.13$coefficients[-1], "p.MLE" = summary(lm.MLE.13)[[4]][-1,4],
                               "est.heur" = lm.heur.13$coefficients[-1], "p.heur" = summary(lm.heur.13)[[4]][-1,4]),3)

R.sq.13 <- data.frame("sim" = summary(lm.sim.13)[[8]], "ind" = summary(lm.ind.13)[[8]], "MLE" = summary(lm.MLE.13)[[8]], "heur" = summary(lm.heur.13)[[8]])


#---- 26 weeks-----------
lm.sim.26 <- lm(tot.correct[[2]]$sim_med ~ ., data= predictors.het)
lm.ind.26 <- lm(tot.correct[[2]]$ind_med ~ ., data = predictors.het)
lm.MLE.26 <- lm(tot.correct[[2]]$MLE ~ ., data = predictors.het)
lm.heur.26 <- lm(tot.correct[[2]]$heur ~., data=predictors.het)

drivers.26 <- round(data.frame("est.sim" = lm.sim.26$coefficients[-1], "p.sim" = summary(lm.sim.26)[[4]][-1,4],
                               "est.ind" = lm.ind.26$coefficients[-1], "p.ind" = summary(lm.ind.26)[[4]][-1,4],
                               "est.MLE" = lm.ind.26$coefficients[-1], "p.MLE" = summary(lm.MLE.26)[[4]][-1,4],
                               "est.heur" = lm.heur.26$coefficients[-1], "p.heur" = summary(lm.heur.26)[[4]][-1,4]),3)

R.sq.26 <- data.frame("sim" = summary(lm.sim.26)[[8]], "ind" = summary(lm.ind.26)[[8]], "MLE" = summary(lm.MLE.26)[[8]], "heur" = summary(lm.heur.26)[[8]])


#----52 weeks -------------
lm.sim.52 <- lm(tot.correct[[3]]$sim_med ~ ., data= predictors.het)
lm.ind.52 <- lm(tot.correct[[3]]$ind_med ~ ., data = predictors.het)
lm.MLE.52 <- lm(tot.correct[[3]]$MLE ~ ., data = predictors.het)
lm.heur.52 <- lm(tot.correct[[3]]$heur ~., data=predictors.het)

drivers.52 <- round(data.frame("est.sim" = lm.sim.52$coefficients[-1], "p.sim" = summary(lm.sim.52)[[4]][-1,4],
                               "est.ind" = lm.ind.52$coefficients[-1], "p.ind" = summary(lm.ind.52)[[4]][-1,4],
                               "est.MLE" = lm.ind.52$coefficients[-1], "p.MLE" = summary(lm.MLE.52)[[4]][-1,4],
                               "est.heur" = lm.heur.52$coefficients[-1], "p.heur" = summary(lm.heur.52)[[4]][-1,4]),3)

R.sq.52 <- data.frame("sim" = summary(lm.sim.52)[[8]], "ind" = summary(lm.ind.52)[[8]], "MLE" = summary(lm.MLE.52)[[8]], "heur" = summary(lm.heur.52)[[8]])


lm.13 <- list("lm.sim.13" = lm.sim.13, "lm.ind.13" = lm.ind.13, "lm.MLE.13" = lm.MLE.13,"lm.heur.13" = lm.heur.13)
lm.26 <- list("lm.sim.26" = lm.sim.26, "lm.ind.26" = lm.ind.26, "lm.MLE.26" = lm.MLE.26,"lm.heur.26" = lm.heur.26)
lm.52 <- list("lm.sim.52" = lm.sim.52, "lm.ind.52" = lm.ind.52, "lm.MLE.52" = lm.MLE.52,"lm.heur.52" = lm.heur.52)

save(lm.13, lm.26, lm.52, file="results/221_active customers regression models.Rdata")  

