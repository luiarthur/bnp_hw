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

BB <- 2000
# DP ###
system.time( sourceCpp("dppois.cpp") )#, showOutput=T) )
system.time(z <- dppois(y,cs(x),a_mu=1,b_mu=1,a_tau=1,b_tau=1,
                        a_alpha=1,b_alpha=1,m_beta=0,s_beta=1,
                        cs_mu=2,cs_beta=.2,B=100000))
z$acc_beta
z$acc_mu

talpha <- tail(z$alpha,BB)
ttheta <- tail(z$theta,BB)
tbeta <- tail(z$beta,BB)
tmu <- tail(z$mu,BB)
ttau <- tail(z$tau,BB)

par(mfrow=c(4,1))
  plot(talpha,main=bquote(alpha),ylab='',type='l')
  plot(tmu,main=bquote(mu),ylab='',type='l')
  plot(ttau,main=bquote(tau),ylab='',type='l')
  plot(tbeta,main=bquote(beta),ylab='',type='l')
par(mfrow=c(1,1))

ut <-  apply(ttheta,1,function(x) length( unique(x) ))
tut <- table(ut)
plot(tut / sum(tut), main=paste("Mean Number of Clusters =",mean(ut)))

plot(ttheta[,1],cex=.1,col=rgb(.1,.1,.1,.2),main="theta")
for ( k in 1:ncol(z$theta) )
  points(ttheta[,k],cex=.1,col=rgb(.1,.1,.1,.2))

plot(tail(z$mu,BB), tail(z$tau,BB),type='p',pch=20)

# Post Pred
y0.pred <- matrix(0,BB,n)
y.pred <- matrix(0,BB,n)
theta.pred <- apply(matrix(1:(BB)),1,function(b) {
            th <- 0

            if (talpha[b] / (talpha[b]+n) > runif(1)) {
              th <- rgamma(1,tmu[b],rate=ttau[b])
            } else {
              th <- sample(ttheta[b,],1)
            }

            th
      })
for (i in 1:n) y0.pred[,i] <- rpois(BB, theta.pred*exp(cs(x)[i]*tbeta) )
for (i in 1:n) y.pred[,i] <- rpois(BB, ttheta[,i]*exp(cs(x)[i]*tbeta) )
ord <- order(x)
plot(x[ord],apply(y0.pred,2,mean)[ord],col="green",pch=20,type='o')
points(x[ord],y[ord],col="grey80",cex=2,pch=20)
points(x[ord],apply(y.pred,2,mean)[ord],col="blue",pch=20)

# Simple Poisson Regression ######
system.time( sourceCpp("regpois.cpp") )
system.time(rgp <- regpois(y,cs(x),zeta=1,mu=1,m=0,s2=.1,cs_beta=.2,B=1000000))
rgp$acc_beta
par(mfrow=c(2,1))
plot(tail(rgp$theta,BB),type='l'); mean(tail(rgp$theta,BB))
plot(tail(rgp$beta,BB),type='l'); mean(tail(rgp$beta,BB))
par(mfrow=c(1,1))

#plot(tail(rgp$beta,BB), tail(rgp$theta,BB), col=rgb(.5,.5,.5,.1:BB/BB),pch=20,cex=.1)
