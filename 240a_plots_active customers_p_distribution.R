load("results/201_active_customers_optimised_p.Rdata")
par(mfrow=c(1,2))
plot(density(opt.p[[1]]$sim.med),col="red", main="optimal p for T = 13 weeks", xlim=c(0,1),
     ylab="",xlab="")
lines(density(opt.p[[1]]$ind.med), col="forestgreen")
lines(density(opt.p[[1]]$MLE), col="blue")
lines(density(opt.p[[1]]$heur), col="grey43")
# plot.new()
legend("topright",col=c("red","forestgreen","blue","grey43"), lty=1, 
       legend=c("sim","ind","MLE","heur"))
