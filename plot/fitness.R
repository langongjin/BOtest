#Gongjin Lan
rm(list = ls())
library(ggplot2)
library(dplyr)

#setwd("~/projects/evert-simulator/tol-revolve/results-rl/spider9/-40/")
#fitness = read.table("fitness.txt")
setwd("~/projects/bayesian/BO/plot/")
fitness = read.table("michalewicz_ma_0.2_gp_0.1.txt")


#---------------------------------------------------------------------------#
#distPro2sum3times =data.frame(distPro2$V1[1:8000] + distPro2$V1[8001:16000] + distPro2$V1[16001:24000]) / 3
fitnessmean <- aggregate(fitness,list(rep(1:(nrow(fitness)%/%20+1),each=20,len=nrow(fitness))),mean)[-1]
fitnessmax <- aggregate(fitness,list(rep(1:(nrow(fitness)%/%20+1),each=20,len=nrow(fitness))),max)[-1]

fitmean = data.frame(gen = rep(1:20),fit = fitnessmean$V1)
fitmax = data.frame(gen = rep(1:20),fit = fitnessmax$V1)

fitmanaverage = aggregate( fit ~ gen, fitmean, mean)
fiteamaxaverage = aggregate( fit ~ gen, fitmax, mean)

bestEA = vector(, 20)
for (i in 1:20)
{
  if(i == 1)
  {
    bestEA[1] = fiteamaxaverage$fit[1]
  }
  else if(bestEA[i-1] < fiteamaxaverage$fit[i])
  {
    bestEA[i] = fiteamaxaverage$fit[i]
  }
  else
    bestEA[i] = bestEA[i-1]
}

ggplot(fitmanaverage, aes(1:nrow(fitmanaverage))) + 
  geom_point(aes(1:nrow(fitmanaverage), fitmanaverage$fit), alpha = 1/2, colour = "black", size = 1) +
  geom_line (aes(1:20, bestEA), alpha = 1/2, colour = "black", size = 1) +
  geom_hline(yintercept=4.687658, linetype="dashed", color = "red", size=0.5) +
  #geom_smooth(aes(1:nrow(fitmanaverage), fitmanaverage$fit), colour = "black", method = "loess") +
  labs(x="Generations", y="fitness",title="michalewicz_ma_0.2_gp_0.1") +
  ylim(0,5) 
#scale_colour_manual("", breaks = c("k=2", "k=5"), values=c("red","black"))+
theme(#legend.position="topleft", 
  #legend.title = element_text("k parameters"),
  #plot.title = element_text(hjust = 0.5),
  legend.text = element_text(colour="blue", size = 16, face = "bold"),
  panel.grid.major = element_line(colour = "grey95", size = 0.25),
  panel.background = element_rect(fill = "transparent", colour="black"))


#-------------------------------------------------------------------------#
#notes: rep(1:500, each=20) mean repeat 20 times for numbers of 1,2,3...500.
#note the difference with rep(1:500,times=20)

fit_for_plot = data.frame(gen = rep(1:20, each=20),fit = fitness$V1)
# get top 3 in same gen by the functions group_by() and top_n()
fit_top3 <- fit_for_plot %>% group_by(gen) %>% top_n(n = 3, wt = fit)

fit_top3 <- unique(fit_top3)

ggplot(fit_for_plot, aes(gen, fit)) + 
  geom_point(alpha = 1/2, colour = "black", size = 0.5) +
  geom_smooth(colour = "blue", method = "loess") +
  geom_smooth(data = fit_top3, colour = "red", method = "loess") +
  labs(x="Generations", y="fitness") +
       #title="The fitness of directed locomotion") +
  ylim(-50000,0) 



