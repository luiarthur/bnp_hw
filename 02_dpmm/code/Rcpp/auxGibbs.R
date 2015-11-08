library(Rcpp)
sourceCpp("dp_arma.cpp")

y <- read.table("../../dat/hw2.dat")[[1]]
#plot(density(y))

B <- 100000
burn <- B * .3
system.time(out <- gibbs(y,1,1,2,burn,B))

