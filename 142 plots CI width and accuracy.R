load("results/102_range_accuracy_ind.Rdata")
range.accuracy.ind <- correct.interquartile.ratio.active.only
range.ind <- relative.range.active
# range.ind.50 <- relative.active.range.75
range.ind.90 <- relative.active.range.90
# range.ind.2sd <- relative.active.range.sd


load("results/103_range_accuracy_mcmc.Rdata")
range.accuracy.mcmc <- correct.interquartile.ratio.active.only
#range.mcmc <- relative.range.active
# range.mcmc.50 <- relative.active.range.75
range.mcmc.90 <- relative.active.range.90
# range.mcmc.2sd <- relative.active.range.sd

# first check
# x11()
# par(mfrow=c(2,2))
# boxplot(range.accuracy.ind, outline=F,main="accuracy ind closed",ylim=c(0,100))
# boxplot(range.accuracy.mcmc, outline=F, main="accuracy mcmc sim",ylim=c(0,100))
# boxplot(range.ind, outline=F, main="rel range ind closed",ylim=c(0,4))
# boxplot(range.mcmc, outline=F, main="rel range mcmc sim",ylim=c(0,4))

load("results/112_ci_accuracy_real.Rdata")
load("dataset_characteristics_real.Rdata")

ci.accuracy.ind.real <- ci.accuracy.ind.real[1:22,]
ci.accuracy.sim.real <- ci.accuracy.sim.real[1:22,]
ci.range.id.real <- ci.range.id.real[1:22,]
ci.range.sim.real <- ci.range.sim.real[1:22,]

x11()
par(mfrow=c(1,2))
boxplot(data.frame("ind"=range.accuracy.ind$X90., "sim" = range.accuracy.mcmc$X90.),outline=F, main="90% CI accuracy")
points(rep(1,13),ci.accuracy.ind.real[dataset.characteristics.real$dataset!=1,2],pch=4,col="red")
points(rep(2,13),ci.accuracy.sim.real[dataset.characteristics.real$dataset!=1,2],pch=4,col="red")

boxplot(data.frame("ind"=range.ind.90, "sim" = range.mcmc.90),main="relative range width",outline=F)
points(rep(1,13),ci.range.id.real[dataset.characteristics.real$dataset!=1,2],pch=4,col="red")
points(rep(2,13),ci.range.sim.real[dataset.characteristics.real$dataset!=1,2],pch=4,col="red")


View(data.frame(dataset.characteristics.real[,1:3],"acc.sim"=round(ci.accuracy.sim.real[,2],1), "acc.ind" = round(ci.accuracy.ind.real[,1],1),"range.sim"=round(ci.range.sim.real[,2],1), "range.ind"=round(ci.range.id.real[,2],1)))

