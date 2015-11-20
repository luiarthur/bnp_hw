# If an older version of gcc is installed globally, but you have a newer 
# version installed locally, add this to ~/.R/Makevars:
# CXX = ~/gcc-5.2.0/bin/g++
# replace with the path to new local g++
library(Rcpp)
source("../../../R_Functions/plotPost.R",chdir=T)
dat <- read.table("../../dat/fabric.dat",header=T)
cs <- function(x) ( x - mean(x) ) / sd(x)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo

# DP ###
system.time( sourceCpp("dppois.cpp") )#, showOutput=T) )
system.time(z <- dppois(y,cs(x),a_mu=1,b_mu=1,a_tau=1,b_tau=1,
                        a_alpha=1,b_alpha=1,m_beta=1,s_beta=1,
                        cs_mu=1,cs_beta=.2,B=100000))
z$acc_beta
z$acc_mu

par(mfrow=c(4,1))
  plot(tail(z$alpha,2000),main=bquote(alpha),ylab='',type='l')
  plot(tail(z$mu,2000),main=bquote(mu),ylab='',type='l')
  plot(tail(z$tau,2000),main=bquote(tau),ylab='',type='l')
  plot(tail(z$beta,2000),main=bquote(beta),ylab='',type='l')
par(mfrow=c(1,1))

ut <-  apply(tail(z$theta,1000),1,function(x) length( unique(x) ))
tut <- table(ut)
plot(tut / sum(tut))

plot(tail(z$theta[,1],2000),cex=.1,col=rgb(.1,.1,.1,.1))
for ( k in 1:ncol(z$theta) )
  points(tail(z$theta[,k],2000),cex=.1,col=rgb(.1,.1,.1,.1))

# Simple Poisson Regression ######
system.time( sourceCpp("regpois.cpp") )
system.time(rgp <- regpois(y,cs(x),a=1,b=1,m=0,s2=.1,cs_beta=.2,B=1000000))
rgp$acc_beta
par(mfrow=c(2,1))
plot(tail(rgp$theta,2000),type='l'); mean(tail(rgp$theta,2000))
plot(tail(rgp$beta,2000),type='l'); mean(tail(rgp$beta,2000))
par(mfrow=c(1,1))

plot(tail(rgp$beta,5000), tail(rgp$theta,5000), col=rgb(.5,.5,.5,.1:5000/5000),pch=20,cex=.1)

