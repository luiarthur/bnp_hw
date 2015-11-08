library(Rcpp)
sourceCpp("dp_arma.cpp")
source("auxGibbs.R")
y <- read.table("../../dat/hw2.dat")[[1]]

B <- 10000
burn <- round(B * .3)

# Rcpp is 10 times faster than R
system.time(outCpp <- auxGibbsCpp(y,a=2,s=1,cs=3,B=B))
system.time(  outR <-   auxGibbsR(y,a=2,s=1,cs=3,B=B))

out <- list(outCpp, outR)

par(mfrow=c(1,2))
for ( i in 1:length(out) ) {
  o <- out[[i]]
  mean.theta <- apply(tail(o,B-burn),2,mean)
  plot(y,pch=20,col='grey',cex=2, main=ifelse(i==1,"Rcpp","R"))
  points(mean.theta,pch=20,col=rgb(0,0,1,.5))
}
par(mfrow=c(1,1))

 plot(density(apply(tail(out[[1]],B-burn),2,mean)), ylim=c(0,.4), lwd=2, col="purple")
lines(density(apply(tail(out[[2]],B-burn),2,mean)), ylim=c(0,.4), lwd=2, col="red")

 plot(density(apply(tail(out[[1]],B-burn),2,sd)), ylim=c(0,3), lwd=2, col="purple")
lines(density(apply(tail(out[[2]],B-burn),2,sd)), ylim=c(0,3), lwd=2, col="red")

