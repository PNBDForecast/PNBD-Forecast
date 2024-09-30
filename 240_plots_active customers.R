load("results/201_active_customers.Rdata")
load("results/212_active_customers_real.Rdata")
load("results/forecasts/forecast_true.Rdata")
load("results/201_active_customers_optimised_p.Rdata")
load("sample_analysis_real.Rdata")

#####################################################
# plot tot.correct and opt.correct
#####################################################
  # t.star definiert den Index für den Plot


# Plot densities of p = 0.5 and optimum threshold  
plot.1 <- function(t.star)  {
    par(mfrow=c(2,2))
    
    plot(density(tot.correct[[t.star]]$sim_med), main="tot.correct sim.med",xlab="",ylab="")
    lines(density(opt.correct[[t.star]]$sim.med), col="red")
    abline(v=mean(tot.correct[[t.star]]$sim_med))
    abline(v=mean(opt.correct[[t.star]]$sim.med),col="red")
    legend("topleft", col=c(1,2),lty=1, legend=c("p = 0.5", "optimised p"))
    
    plot(density(tot.correct[[t.star]]$ind_med), main="tot.correct ind median",xlab="",ylab="")
    lines(density(opt.correct[[t.star]]$ind.med), col="red")
    abline(v=mean(tot.correct[[t.star]]$ind_med))
    abline(v=mean(opt.correct[[t.star]]$ind.med),col="red")
    
    plot(density(tot.correct[[t.star]]$MLE, na.rm=T), main="tot.correct MLE",xlab="",ylab="",ylim=c(0,0.08))
    lines(density(opt.correct[[t.star]]$MLE, na.rm=T), col="red")
    abline(v=mean(tot.correct[[t.star]]$MLE, na.rm=T))
    abline(v=mean(opt.correct[[t.star]]$MLE, na.rm=T),col="red")
    
    plot(density(tot.correct[[t.star]]$heur), main="tot.correct heur",xlab="",ylab="", ylim=c(0,0.08))
    lines(density(opt.correct[[t.star]]$heur), col="red")
    abline(v=mean(tot.correct[[t.star]]$heur, na.rm=T))
    abline(v=mean(opt.correct[[t.star]]$heur, na.rm=T),col="red")
    }
  
  
# Scatterplot
plot.2 <- function(){
    par(mfrow=c(1,1))
    plot(sensitivity[[3]]$heur, precision[[3]]$heur, col="gray87",ylim=c(0,100), xlim=c(0,100), main="all data sets",xlab="",ylab="",pch=16)
    points(sensitivity[[3]]$MLE, precision[[3]]$MLE, col="lightblue",pch=16)
    points(sensitivity[[3]]$sim_med, precision[[3]]$sim_med, col="lightpink",pch=16)
    
    points(sensitivity.real[[3]]$heur[sample.analysis.real$dataset !=1], precision.real[[3]]$heur[sample.analysis.real$dataset !=1], col="gray31",pch=16)
    points(sensitivity.real[[3]]$MLE[sample.analysis.real$dataset !=1], precision.real[[3]]$MLE[sample.analysis.real$dataset !=1], col="blue",pch=16)
    points(sensitivity.real[[3]]$sim_med[sample.analysis.real$dataset !=1], precision.real[[3]]$sim_med[sample.analysis.real$dataset !=1], col="red",pch=16)
    
    
    legend("bottomleft",legend=c("ind.median","MLE","heur"),col=c("red","blue","gray31"),bty="n",pch=rep(16,3))
  }
  

  t.star.value <- c(13,26,52)
  plot.3 <- function(t.star){
    boxplot(tot.correct[[t.star]], main=paste("correctly assigned customers T* = ",t.star.value[t.star],"weeks"), outline=FALSE)
    points(rep(1,13), tot.correct.real[[t.star]]$sim_med[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(2,13), tot.correct.real[[t.star]]$ind_mean[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(3,13), tot.correct.real[[t.star]]$het_med[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(4,13), tot.correct.real[[t.star]]$sim_mean[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(5,13), tot.correct.real[[t.star]]$ind_mean[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(6,13), tot.correct.real[[t.star]]$het_mean[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(7,13), tot.correct.real[[t.star]]$MLE[sample.analysis.real$dataset !=1], col="red",pch = 4)
    points(rep(8,13), tot.correct.real[[t.star]]$heur[sample.analysis.real$dataset !=1], col="red",pch = 4)
    
  }
  
 
# Vergleich aller simulierten Datensätze mit denjenigen, bei denen der Anteil der active buyers > 10% ist.  
  
  
  # Typ 1: Boxplots
  active.customer.ratio <- data.frame("13" = rep(NA,3000), "26"=rep(NA,3000), "52" = rep(NA,3000))
  for (row in 1:3000){
    for (tstar in 1:3){
      active.customer.ratio[row,tstar] <- 100*(sum(x.star.true[[row]][,tstar] > 0) / nrow(x.star.true[[row]]))
    }
  }
  
  # change missing values from 0 to NA
  tot.correct[[3]]$MLE[tot.correct[[3]]$MLE == 0] <- NA
  tot.correct[[2]]$MLE[tot.correct[[2]]$MLE == 0] <- NA
  tot.correct[[1]]$MLE[tot.correct[[1]]$MLE == 0] <- NA
  
  
  par(mfrow=c(3,2))
  boxplot(tot.correct[[3]]$sim_med[active.customer.ratio[,3] > 10], main="sim all data sets",ylim=c(70,100))
  boxplot(tot.correct[[3]]$sim_med, main="sim > 10% active buyers only",ylim=c(70,100))
  boxplot(tot.correct[[3]]$MLE[active.customer.ratio[,3] > 10], main="sim all data sets",ylim=c(50,100))
  boxplot(tot.correct[[3]]$MLE, main="sim > 10% active buyers only",ylim=c(50,100))
  boxplot(tot.correct[[3]]$heur[active.customer.ratio[,3] > 10], main="sim all data sets",ylim=c(40,100))
  boxplot(tot.correct[[3]]$heur, main="sim > 10% active buyers only",ylim=c(40,100))  
  
  
  # Typ 2: Scatter Plots
  par(mfrow=(c(1,2)))
  plot(sensitivity[[3]]$heur, precision[[3]]$heur, col="gray87",ylim=c(0,100), xlim=c(0,100), main="all data sets",xlab="",ylab="",pch=16)
  points(sensitivity[[3]]$MLE, precision[[3]]$MLE, col="lightblue",pch=16)
  points(sensitivity[[3]]$sim_med, precision[[3]]$sim_med, col="lightpink",pch=16)
  
  plot(sensitivity[[3]]$heur[active.customer.ratio[,3] > 10], precision[[3]]$heur[active.customer.ratio[,3] > 10], col="gray87",ylim=c(0,100), xlim=c(0,100), main="data sets > 10% active",xlab="",ylab="",pch=16)
  points(sensitivity[[3]]$MLE[active.customer.ratio[,3] > 10], precision[[3]]$MLE[active.customer.ratio[,3] > 10], col="lightblue",pch=16)
  points(sensitivity[[3]]$sim_med[active.customer.ratio[,3] > 10], precision[[3]]$sim_med[active.customer.ratio[,3] > 10], col="lightpink",pch=16)
  