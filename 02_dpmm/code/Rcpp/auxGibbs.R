library(Rcpp)
sourceCpp("dp_arma.cpp")

y <- read.table("../../dat/hw2.dat")[[1]]
#plot(density(y))

B <- 1000
burn <- B * .3

system.time(out <- gibbs(y,1,1,2,B)) # 30 times faster than julia

tab <- table(out[B,])
mt <- apply(out,2,mean)
mns <- apply(tail(out,B-burn),1,function(x) length(unique(x)))
tabm <- table(mns) / sum(table(mns))
plot(tabm)
