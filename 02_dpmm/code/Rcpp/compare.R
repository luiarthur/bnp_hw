library(Rcpp)
sourceCpp("dp_arma.cpp")
source("auxGibbs.R")
y <- read.table("../../dat/hw2.dat")[[1]]

B <- 10000
burn <- B * .3


system.time(out <- gibbs(y,a=2,s=1,cs=3,B=10000)) # 30 times faster than julia

tab <- table(out[B,])
mt <- apply(tail(out,B-burn),2,mean)
mns <- apply(tail(out,B-burn),1,function(x) length(unique(x)))
tabm <- table(mns) / sum(table(mns))
plot(tabm)

plot(y,pch=20,col='grey',cex=2)
points(mt,pch=20,col=rgb(0,0,1,.5))
