load("results/forecasts_topA.Rdata")

table(fc.true[[1]])
Q90 <- unlist(lapply(fc.true, function(x) quantile(x,0.9, na.rm=TRUE)))
Q80 <- unlist(lapply(fc.true, function(x) quantile(x,0.8, na.rm=TRUE)))


active.cust <- unlist(lapply(fc.true, function(x) sum(x > 0) / length(x)))

many.q90 <- unlist(lapply(fc.true, function(x) sum(x >= quantile(x,0.9, na.rm=TRUE))/length(x)))
many.q80 <- unlist(lapply(fc.true, function(x) sum(x >= quantile(x,0.8, na.rm=TRUE))/length(x)))


table(many.q90[active.cust > 0.1])
table(many.q80[active.cust > 0.1])


# if we remove those data sets where the ratio of active customers is < 10%, we have 1,716 data sets left where
# the largest ratio of customers with repurchases 
# 
Top.10.accuracy <- matrix(nrow=3000, ncol= 4) # accuracy summary
colnames(Top.10.accuracy) <- c("MCMC median","MCMC mean", "MLE", "heur")
top.10.count <- Top.20.accuracy <- Top.10.accuracy <- as.data.frame(Top.10.accuracy)
top.10.count$true <- NA
top.20.count <- top.10.count


for (row in 1:3000){
  
  
  #-----------------------------------------
  top.10.true <- which (fc.true[[row]] >= quantile(fc.true[[row]],0.9, na.rm=TRUE) ) # may be more than 10%!
  top.10.heur <- order (fc.heur[[row]], decreasing =TRUE)[1:length(top.10.true)]
  top.10.ind <- order (fc.ind[[row]], decreasing =TRUE)[1:length(top.10.true)]
  top.10.ind.mean <- order (fc.ind.mean[[row]], decreasing =TRUE)[1:length(top.10.true)]
  if (!is.null(fc.MLE[[row]])) top.10.MLE <- order (fc.MLE[[row]], decreasing =TRUE)[1:length(top.10.true)]
  
  top.10.count$`MCMC median`[row] <- length(top.10.ind) / (length(fc.true[[row]]))
  top.10.count$`MCMC mean`[row] <- length(top.10.ind.mean) / (length(fc.true[[row]]))
  top.10.count$MLE[row] <- length(top.10.MLE) / (length(fc.true[[row]]))
  top.10.count$heur[row] <- length(top.10.heur) / (length(fc.true[[row]]))
  top.10.count$true[row] <- length(top.10.true) / (length(fc.true[[row]]))
  
  top.20.true <- which (fc.true[[row]] >= quantile(fc.true[[row]],0.95, na.rm=TRUE) )
  top.20.heur <- order (fc.heur[[row]], decreasing =TRUE)[1:length(top.20.true)]
  top.20.ind <- order (fc.ind[[row]], decreasing =TRUE)[1:length(top.20.true)]
  top.20.ind.mean <- order (fc.ind.mean[[row]], decreasing =TRUE)[1:length(top.20.true)]
  if (!is.null(fc.MLE[[row]])) top.20.MLE <- order (fc.MLE[[row]], decreasing =TRUE)[1:length(top.20.true)]
  
  top.20.count$`MCMC median`[row] <- length(top.20.ind) / (length(fc.true[[row]]))
  top.20.count$`MCMC mean`[row] <- length(top.20.ind.mean) / (length(fc.true[[row]]))
  if (!is.null(fc.MLE[[row]])) top.20.count$MLE[row] <- length(top.20.MLE) / (length(fc.true[[row]]))
  top.20.count$heur[row] <- length(top.20.heur) / (length(fc.true[[row]]))
  top.20.count$true[row] <- length(top.20.true) / (length(fc.true[[row]]))
  
  #determine top A% accuracy
  Top.10.accuracy$`MCMC median`[row] <- length(intersect(top.10.true, top.10.ind)) / (length(top.10.true))
  Top.10.accuracy$`MCMC mean`[row] <-  length(intersect(top.10.true, top.10.ind.mean)) / (length(top.10.true))
  Top.10.accuracy$MLE[row] <-  length(intersect(top.10.true, top.10.MLE)) / (length(top.10.true))
  Top.10.accuracy$heur[row] <-  length(intersect(top.10.true, top.10.heur)) / (length(top.10.true))
  
  Top.20.accuracy$`MCMC median`[row] <-  length(intersect(top.20.true, top.20.ind)) / (length(top.20.true))
  Top.20.accuracy$`MCMC mean`[row] <-  length(intersect(top.20.true, top.20.ind.mean))/ (length(top.20.true))
  if (!is.null(fc.MLE[[row]])) Top.20.accuracy$MLE[row] <-  length(intersect(top.20.true, top.20.MLE)) / (length(top.20.true))
  Top.20.accuracy$heur[row] <-  length(intersect(top.20.true, top.20.heur))/ (length(top.20.true))
  
  
  
}
save(Top.10.accuracy, Top.20.accuracy, file="results/TopA_accuracy.Rdata")

boxplot.data <- as.data.frame(rbind(cbind("MCMC",Top.10.accuracy$`MCMC median`), cbind("MLE",Top.10.accuracy$MLE),cbind("heuristic",Top.10.accuracy$heur)))
colnames(boxplot.data) <- c("method","Top 10% accuracy")
boxplot.data$"Top 20% accuracy" <- c(Top.20.accuracy$`MCMC median`, Top.20.accuracy$MLE, Top.20.accuracy$heur)
boxplot.data$`Top 10% accuracy` <- as.numeric(boxplot.data$`Top 10% accuracy`)
boxplot.data$`Top 20% accuracy` <- as.numeric(boxplot.data$`Top 20% accuracy`)

boxplot.data$method <- factor(boxplot.data$method, levels=c("MCMC","MLE","heuristic"))
x11()
par(mfrow=c(2,1))
boxplot(boxplot.data$`Top 10% accuracy` ~ boxplot.data$method, col=c("red","blue","gray31"),ylab="",xlab="",main="Top 10% accuracy")
boxplot(boxplot.data$`Top 20% accuracy` ~ boxplot.data$method, col=c("red","blue","gray31"),ylab="",xlab="",main="Top 20% accuracy")

