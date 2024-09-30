
##############################################################
# individual deviation prediction
##############################################################
# Die prediction hat, unabhängig vom Verfahren (mcmc, ind, MLE), der forecast period und der 
# Art der Predictor (deskriptiv oder Parameterschätzer) ein R² < 0.4 und ist daher unbrauchbar.

# Dies gilt ebenfalls bei der Einschränkung auf repeat buyer (sowohl die absolute deviation als auch die relative).


# bereits gerechnet und gespeichert: x.stars als Vektoren ("results/x_star_vectors.Rdata")


create.x.vectors <- function(){
  load("dataset_characteristics.Rdata")
  load("results/forecast_MLE.Rdata")
  load("results/forecast_true.Rdata")
  load("results/forecast_mcmc.Rdata")
  
  x.star.true.13 <- x.star.mcmc.13 <- x.star.ind.13 <- x.star.MLE.13 <- 0
  x.star.true.26 <- x.star.mcmc.26 <- x.star.ind.26 <- x.star.MLE.26 <- 0
  x.star.true.52 <- x.star.mcmc.52 <- x.star.ind.52 <- x.star.MLE.52 <- 0
  
  for (row in 1:3000){
    x.star.true.13 <- c(x.star.true.13, x.star.true[[row]][,1])
    x.star.true.26 <- c(x.star.true.26, x.star.true[[row]][,2])
    x.star.true.52 <- c(x.star.true.52, x.star.true[[row]][,3])
    
    x.star.mcmc.13 <- c(x.star.mcmc.13, x.star.mcmc[[row]][1,,1,2])
    x.star.mcmc.26 <- c(x.star.mcmc.26, x.star.mcmc[[row]][2,,1,2])
    x.star.mcmc.52 <- c(x.star.mcmc.52, x.star.mcmc[[row]][3,,1,2])
    
    x.star.ind.13 <- c(x.star.ind.13, x.star.mcmc[[row]][1,,2,2])
    x.star.ind.26 <- c(x.star.ind.26, x.star.mcmc[[row]][2,,2,2])
    x.star.ind.52 <- c(x.star.ind.52, x.star.mcmc[[row]][3,,2,2])
    
    x.star.MLE.13 <- c(x.star.MLE.13, x.star.MLE[[row]][,1])
    x.star.MLE.26 <- c(x.star.MLE.26, x.star.MLE[[row]][,2])
    x.star.MLE.52 <- c(x.star.MLE.52, x.star.MLE[[row]][,3])
    
    cat(format(Sys.time(), '%H:%M:%S'),"Datensatz",row,"beendet.","\n")
  }
  
  x.star.mcmc.13 <- x.star.mcmc.13[-1]
  x.star.mcmc.26 <- x.star.mcmc.26[-1]
  x.star.mcmc.52 <- x.star.mcmc.52[-1]
  x.star.true.13 <- x.star.true.13[-1]
  x.star.true.26 <- x.star.true.26[-1]
  x.star.true.52 <- x.star.true.52[-1]
  x.star.ind.13 <- x.star.ind.13[-1]
  x.star.ind.26 <- x.star.ind.26[-1]
  x.star.ind.52 <- x.star.ind.52[-1]
  x.star.MLE.13 <- x.star.MLE.13[-1]
  x.star.MLE.26 <- x.star.MLE.26[-1]
  x.star.MLE.52 <- x.star.MLE.52[-1]
  
  rm(x.star.mcmc, x.star.true, x.star.MLE)
  
  save.image("helpers/x_star_vectors.Rdata")
}
