set.seed(1105)
sample.analysis <- data.frame("cohort.size" = round(runif(3000,1000,4000)),
                              "calibration" = round(runif(3000,26,72)),
                              "E(lambda)" = runif(3000, 0.02, 0.3),
                              "CV(lambda)" = runif(3000, 0.5, 2.5),
                              "E(mu)" = runif(3000, 0.02, 0.2),
                              "CV(mu)" = runif(3000, 0.5, 2.5))

sample.analysis$r <- sample.analysis$CV.lambda.^(-2)
sample.analysis$alpha <- sample.analysis$r / sample.analysis$E.lambda.
sample.analysis$s <- sample.analysis$CV.mu.^(-2)
sample.analysis$beta <- sample.analysis$s / sample.analysis$E.mu. 

save(sample.analysis, file="sample_analysis.Rdata")