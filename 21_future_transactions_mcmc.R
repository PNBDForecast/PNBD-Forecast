# draws future transaction based on MCMC parameter draws

# wird keine sample_size vorgegeben, dann wird die sample_size der Parameterschätzung übernommen

library(data.table)
library(parallel)
library(BTYD)
library(coda)
library(mcmcplots)
library(plyr)


draw_future_transactions <- function(cal.cbs,draws,T.star=cal.cbs$T.star,sample_size=NULL){
  
  
  # Festlegung der sample size
  if (is.null(sample_size)) {
    nr_of_draws <- niter(draws$level_2) * nchain(draws$level_2)
  } else {
    stopifnot(is.numeric(sample_size))
    nr_of_draws <- as.integer(sample_size)
  }
  
  # Festlegung der Dimensionen und Parameternamen
  nr_of_cust <- dim(draws)[3]
  parameters <- varnames(dim(draws[3]))
  
  
  
  # Erstellung der ystar-Matrix
  x.stars <- array(NA_real_, dim = c(nr_of_draws, nr_of_cust))
  
  
  
  # Formatierung von T.star (falls erforderlich)
  if (length(T.star) == 1)
    T.star <- rep(T.star, nr_of_cust)
  
  
  # Definition der truncated exponential distribution
  
  draw_left_truncated_exp <- function(lower, lambda) {
    rand <- runif(1, pexp(lower, lambda), 1)
    qexp(rand, lambda)
  }
  
  ## Schleife für die Simulation
  
  for (cust in 1:nrow(cal.cbs)) {
    
    # Definition der festen Parameter Tcal, Tstar, tx
    Tcal <- cal.cbs[cust,3]
    Tstar <- T.star[cust]
    tx <- cal.cbs[cust,2]
    
    
    # Einlesen der simulierten Parameter tau und lambda
    taus <- MCMC_ind[3,,cust]
    lambdas <- MCMC_ind[1,,cust]
    stopifnot(length(taus) == length(lambdas))
    if (!is.null(sample_size)) {
      idx <- sample(length(taus), size = sample_size, replace = TRUE)
      taus <- taus[idx]
      lambdas <- lambdas[idx]
    }
    alive <- (taus > Tcal)
    
    # Case: customer alive
    for (draw in which(alive)) {
      # sample itt which is larger than (Tcal-tx)
      itts <- draw_left_truncated_exp(Tcal - tx, lambdas[draw])
      # sample 'sufficiently' large amount of inter-transaction times
      minT <- pmin(Tcal + Tstar - tx, taus[draw] - tx)
      nr_of_itt_draws <- pmax(10, round(minT * lambdas[draw]))
      itts <- c(itts, rexp(nr_of_itt_draws * 2, lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rexp(nr_of_itt_draws * 4, rate = lambdas[draw]))
      if (sum(itts) < minT)
        itts <- c(itts, rexp(nr_of_itt_draws * 800, rate = lambdas[draw]))
      if (sum(itts) < minT)
        stop("not enough inter-transaction times sampled! cust:", cust, " draw:", draw, " ", sum(itts),
             " < ", minT)
      x.stars[draw, cust] <- sum(cumsum(itts) < minT)
    }
    
    # Case: customer churned
    if (any(!alive)) {
      x.stars[!alive, cust] <- 0
    }
  }
  return(x.stars)

}




