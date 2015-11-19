library(Rcpp)
library(RcppArmadillo) # For making packages
source("../../../R_Functions/plotPost.R",chdir=T)
dat <- read.table("../../dat/fabric.dat",header=T)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11
#Sys.setenv("PKG_CXXFLAGS"="-v") # check version
system.time( sourceCpp("dppois.cpp") )#, showOutput=T) )
system.time(z <- dppois(y,x,a_mu=1,b_mu=1,a_tau=1,b_tau=1,
                        a_alpha=1,b_alpha=1,m_beta=1,s_beta=1,
                        cs_mu=1,cs_beta=.0003,B=50000))

par(mfrow=c(4,1))
plot.post(tail(z$alpha,2000),main=bquote(alpha),ylab='')
plot.post(tail(z$mu,2000),main=bquote(mu),ylab='')
plot.post(tail(z$tau,2000),main=bquote(tau),ylab='')
plot.post(tail(z$beta,2000),main=bquote(beta),ylab='')
par(mfrow=c(4,1))

z$acc_beta
z$acc_mu


