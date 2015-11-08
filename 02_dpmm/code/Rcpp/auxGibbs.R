library(Rcpp)
sourceCpp("dp_arma.cpp")

y <- read.table("../../dat/hw2.dat")[[1]]
#plot(density(y))

B <- 10000
burn <- B * .3

system.time(out <- gibbs(y,1,1,2,B)) # 30 times faster than julia

tab <- table(out[B,])
mt <- apply(out,2,mean)
mns <- apply(out[tail(B-burn),],2,function(x) length(unique(x)))
tabm = table(mns)
sort_table(tabm,true)

