setwd('/Users/avinashbarnwal/Desktop/Research/PhD/network_elasticnet/NYSDS-2019/result/')
library(ggplot2)

data = read.table("AIC_performance.csv",header=TRUE,sep=",")
p    = ggplot(data=data,aes(x=lambda,y=aic,col=factor(alpha))) +
       geom_line() + xlab("Lambda") + ylab("AIC") + theme(legend.position=c(0.8,0.3)) + labs(col="Alpha")

p
