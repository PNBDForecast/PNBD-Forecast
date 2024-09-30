sourceCpp("slice_sampler.cpp")

#============================================================================================================= 

# define log-likelihood function
ll <- function(pa) {
  r <- pa[1]
  a <- pa[2]
  s <- pa[3]
  b <- pa[4]
  
  
  x <- cbs[,1]
  tx <- cbs[,2]
  T.cal <- cbs[,3]
  
  
  # Define the LogLikelihood for the Pareto/NBD
  if(a > b){
    A1 <- ((a+T.cal)/(b+T.cal))^s + s/(r+s+x) * (hyperg_2F1(r+s+x, s+1, r+s+x+1, (a-b)/(a+tx)) * ((a+T.cal)/(a+tx))^(r+s+x) - hyperg_2F1(r+s+x, s+1, r+s+x+1, (a-b)/(a+T.cal)))
    
    LL <- sum(lgamma(r+x) + r*log(a) + s*log(b) - lgamma(r) + log(A1) - (r+x+s)*log(a+T.cal))
    
  } else if(a < b){
    A1 <- ((b+T.cal)/(a+T.cal))^(r+x) + s/(r+s+x) * (hyperg_2F1(r+s+x, r+x, r+s+x+1, (b-a)/(b+tx)) * ((b+T.cal)/(b+tx))^(r+s+x) - hyperg_2F1(r+s+x, r+x, r+s+x+1, (b-a)/(b+T.cal)))
    LL <- sum(lgamma(r+x) + r*log(a) + s*log(b) - lgamma(r) + log(A1)  - (r+s+x) * log(b+T.cal)
    )
  } else{ # a==b
    LL <- sum(
      (lgamma(r+x) + r*log(a) + s*log(b) - lgamma(r)) +
        log(
          (1 / ((a+T.cal)^(r+x) * (b+T.cal)^s)) +
            (s / (r+s+x)) *
            (1 / (a+tx)^(r+s+x) -
               1 / (b+T.cal)^(r+s+x))
        )
    )
    
  }
  return(LL)
}



algorithm <- function(cbs,hyper_prior){
  
#--------
  
  draw_lambda_pnbd_DA <- function(cbs, level_1, level_2) {
    N <- nrow(cbs)
    x <- cbs[,1]
    T.cal <- cbs[,3]
    tau <- level_1["tau", ]
    r <- level_2["r"]
    alpha <- level_2["alpha"]
    
    lambda <- rgamma(n = N, shape = r + x, rate = alpha + pmin(tau, T.cal))
    lambda[lambda == 0 | log(lambda) < -30] <- exp(-30)  # avoid numeric overflow
    return(lambda)
  }
  
#--------  
  
  draw_mu_pnbd_DA <- function(cbs, level_1, level_2) {
    N <- nrow(cbs)
    mu <- numeric(N)
    tau <- level_1["tau", ]
    z <- level_1["z", ]
    s <- level_2["s"]
    beta <- level_2["beta"]
    Tcal <- cbs[,3]
    
    mu <- rgamma(n = N, shape = s + 1, rate = beta + tau)
    mu[mu == 0 | log(mu) < -30] <- exp(-30)  # avoid numeric overflow
    return(mu)
    }
  
#--------  
  
  draw_tau <- function(cbs, level_1) {
    N <- nrow(cbs)
    tx <- cbs[,2]
    Tcal <- cbs[,3]
    lambda <- level_1["lambda", ]
    mu <- level_1["mu", ]
    
    mu_lam <- mu + lambda
    t_diff <- Tcal - tx
    
    #draw z 
    p_alive <- 1 / (1 + (mu / mu_lam) * (exp(mu_lam * t_diff) - 1))
    alive <- p_alive > runif(n = N)
    
    #draw tau
    tau <- numeric(N)
    
    # Case: still alive - left truncated exponential distribution -> [Tcal, Inf]
    if (any(alive)) {
      tau[alive] <- Tcal[alive] + rexp(sum(alive), mu[alive])
    }
    
    # Case: churned - double truncated exponential distribution -> [tx, Tcal]
    if (any(!alive)) {
      mu_lam_tx <- pmin(700, mu_lam[!alive] * tx[!alive])
      mu_lam_Tcal <- pmin(700, mu_lam[!alive] * Tcal[!alive])
      # sample with http://en.wikipedia.org/wiki/Inverse_transform_sampling
      rand <- runif(n = sum(!alive))
      tau[!alive] <- -log((1 - rand) * exp(-mu_lam_tx) + rand * exp(-mu_lam_Tcal)) / mu_lam[!alive] # nolint
    }
    
    return(tau)
  }
  
#--------  
  
  
  draw_gamma_params <- function(type, level_1, level_2, hyper_prior) {
    
    #r and alpha
    if (type == "lambda") {
      x <- level_1["lambda", ]
      cur_params <- c(level_2["r"], level_2["alpha"])
      hyper <- as.vector(hyper_prior[1:4])
    } else if (type == "mu") {
      
      #Draw s and beta
      x <- level_1["mu", ]
      cur_params <- c(level_2["s"], level_2["beta"])
      hyper <- as.vector(hyper_prior[5:8])
    }
    slice_sampler_gamma_parameters(x, cur_params, hyper, steps = 5, w = 0.1)
  }
  
#=============================================================================================================

  
  pnbd_DrawParameters <- function(cal.cbs, param_init = NULL, trace = 1000)
  {
    #initialize arrays
    nr_of_cust <- nrow(cbs)
    nr_of_draws <- (mcmc - 1) %/% thin + 1
    level_2_draws <- array(NA_real_, dim = c(nr_of_draws, 5))
    dimnames(level_2_draws)[[2]] <- c("r", "alpha", "s", "beta","LL")
    level_1_draws <- array(NA_real_, dim = c(nr_of_draws, 4, nr_of_cust))
    dimnames(level_1_draws)[[2]] <- c("lambda", "mu", "tau", "z")
    
    #initialize parameters
    level_2 <- level_2_draws[1,1:4]
    level_2["r"] <- runif(1)
    level_2["alpha"] <- runif(1,5,20)
    level_2["s"] <- runif(1)
    level_2["beta"] <- runif(1,5,20)
    #----------
    level_1 <- level_1_draws[1, , ] # nolint
    level_1["lambda", ] <- mean(cbs[,1]) / mean(ifelse(cbs[,2] == 0, cbs[,3], cbs[,2]))
    level_1["tau", ] <- cbs[,2] + 0.5 / level_1["lambda", ]
    level_1["z", ] <- as.numeric(level_1["tau", ] > cbs[,3])
    level_1["mu", ] <- 1 / level_1["tau", ]
    

    run_single_chain <- function(chain_id = 1, cbs, hyper_prior) {
      
      sourceCpp("slice_sampler.cpp")
      set.seed(chain_id)
      
      for (step in 1:(burnin + mcmc)) {
        
        if ( (step - burnin) > 0 & (step - 1 - burnin) %% thin == 0){
                idx <- (step - 1 - burnin) %/% thin + 1
                level_1_draws[idx, , ] <- level_1 # nolint
                level_2_draws[idx, ] <- c(level_2,ll(level_2))
              }
        
        level_1["lambda", ] <- draw_lambda_pnbd_DA(cbs, level_1, level_2)
        level_1["mu", ] <- draw_mu_pnbd_DA(cbs, level_1, level_2)
        level_1["tau", ] <- draw_tau(cbs, level_1)
        level_1["z", ] <- as.numeric(level_1["tau", ] > cbs[,3])
        
        level_2[c("r", "alpha")] <- draw_gamma_params("lambda", level_1, level_2, hyper_prior)
        level_2[c("s", "beta")] <- draw_gamma_params("mu", level_1, level_2, hyper_prior)
      } # end of mcmc steps
      
      return(list(
        "level_1" = lapply(1:nr_of_cust,
                           function(i) mcmc(level_1_draws[, , i], start = burnin, thin = thin)), # nolint
        "level_2" = mcmc(level_2_draws, start = burnin, thin = thin)))
    } # end of run_single_chain
    
    
    #parallel call

    draws <- mclapply(1:chains, function(i) run_single_chain(i, cal.cbs, hyper_prior), mc.cores = ncores)
    
    
    # merge chains into code::mcmc.list objects
    out <- list(level_1 = lapply(1:nrow(cal.cbs), function(i) mcmc.list(lapply(draws, function(draw) draw$level_1[[i]]))),
                level_2 = mcmc.list(lapply(draws, function(draw) draw$level_2)))
    
    return(out)
  } #end of pnbd_DrawParameters
  
  #calling pnbd_DrawParameters
  draws_pnbd <- pnbd_DrawParameters(cbs)
  
  
  return(draws_pnbd)
} # end of algorithm


################
alldraws <- algorithm(cbs,hyper_prior)
