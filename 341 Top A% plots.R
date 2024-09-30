load("results/301_Top_A%.Rdata")
load("results/311_Top_A%_real.Rdata")  

real.idx <- sample.analysis.real$dataset != 1

top.10.df <- data.frame(correct.top.10.ind, correct.top.10.MLE, correct.top.10.heur)
top.20.df <- data.frame(correct.top.20.ind, correct.top.20.MLE, correct.top.20.heur)

par(mfrow=c(2,1))
boxplot(top.10.df[active.customer.ratio > 10,], main="correctly identified Top 10%",outline=F)
points(rep(1,sum(real.idx)), correct.top.10.ind.real[real.idx], col="red",pch=4)
points(rep(2,sum(real.idx)), correct.top.10.MLE.real[real.idx], col="red",pch=4)
points(rep(3,sum(real.idx)), correct.top.10.heur.real[real.idx], col="red",pch=4)


boxplot(top.20.df[active.customer.ratio > 20,], main = "correctly identified Top 20%", outline=F)
points(rep(1,sum(real.idx)), correct.top.20.ind.real[real.idx], col="red",pch=4)
points(rep(2,sum(real.idx)), correct.top.20.MLE.real[real.idx], col="red",pch=4)
points(rep(3,sum(real.idx)), correct.top.20.heur.real[real.idx], col="red",pch=4)
