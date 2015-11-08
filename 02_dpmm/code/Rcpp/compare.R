library(microbenchmark)
library(Rcpp)
sourceCpp("dp_arma.cpp")
source("auxGibbs.R")
y <- read.table("../../dat/hw2.dat")[[1]]

B <- 10000
burn <- B * .3

# Rcpp is 30 times faster than julia
microbenchmark(outCpp <- auxGibbsCpp(y,a=2,s=1,cs=3,B=B),
               outR <- auxGibbsR(y,a=2,s=1,cs=3,B=B))

out <- list(outCpp, outR)[[1]] # change this
mean.theta <- apply(tail(out,B-burn),2,mean)
plot(y,pch=20,col='grey',cex=2)
points(mean.theta,pch=20,col=rgb(0,0,1,.5))
