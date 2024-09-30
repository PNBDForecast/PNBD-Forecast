load("results/120_x_star_regression_models.Rdata")
"sim" <- summary(lm.52$lm.sim.52)[[4]][-1,c(1,4)]
"ind" <- summary(lm.52$lm.ind.52)[[4]][-1,c(1,4)]
"MLE" <- summary(lm.52$lm.MLE.52)[[4]][-1,c(1,4)]
"heur" <- summary(lm.52$lm.heur.52)[[4]][-1,c(1,4)]

coef.table <- round(data.frame("est.sim" = sim[,1], "p.sim" = sim[,2],
                         "est.ind" = ind[,1], "p.ind" = ind[,2],
                         "est.MLE" = MLE[,1], "p.MLE" = MLE[,2],
                         "est.heur" = heur[,1], "p.heur" = heur[,2]),3)
