library(microbenchmark)
library(Rcpp)
sourceCpp("dp_arma.cpp")
source("auxGibbs.R")
y <- read.table("../../dat/hw2.dat")[[1]]

B <- 5000
burn <- B * .3

# Rcpp is 10 times faster than R
microbenchmark(outCpp <- auxGibbsCpp(y,a=2,s=1,cs=3,B=B),
                 outR <-   auxGibbsR(y,a=2,s=1,cs=3,B=B), 
                 times=1, control=list("order"="inorder"))

out <- list(outCpp, outR)

par(mfrow=c(1,2))
for ( i in 1:length(out) ) {
  o <- out[[i]]
  mean.theta <- apply(tail(o,B-burn),2,mean)
  plot(y,pch=20,col='grey',cex=2, main=ifelse(i==1,"Rcpp","R"))
  points(mean.theta,pch=20,col=rgb(0,0,1,.5))
}
par(mfrow=c(1,1))

plot(density(apply(out[[1]],2,mean)),  ylim=c(0,.4), col="purple", lwd=2)
lines(density(apply(out[[2]],2,mean)), ylim=c(0,.4), col="red", lwd=2)

plot(density(apply(out[[1]],2,sd)),  ylim=c(0,3), col="purple", lwd=2)
lines(density(apply(out[[2]],2,sd)), ylim=c(0,3), col="red", lwd=2)

