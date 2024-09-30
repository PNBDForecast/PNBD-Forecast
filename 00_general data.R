library(abind)
library(BTYD)
library(coda)
library(data.table)
library(gsl)
library(mcmcplots)
require(NlcOptim)
library(openxlsx)
library(parallelsugar)
library(plyr)
library(Rcpp)


burnin        <- 2000
ncores        <- 4
mcmc          <- 10000/ncores
thin          <- 20
chains        <- ncores


##############################################
# DEFINE INDIVIDUAL FORECAST IN CLOSED FORM
##############################################
EYt_ind <- function(lambda,mu,z,T.star){
  estimate <- rep(0,times=length(lambda))
  estimate[z==1]  <- (lambda/mu*(1-exp(-mu*T.star)))[z==1]
  return(estimate)
}


# hyper priors
CV_value <- 0.5 # according to Worth the effort study 1
