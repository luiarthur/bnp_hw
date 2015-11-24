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
system.time(z <- dppois(y,cs(x),a_zeta=1,b_zeta=1,a_mu=1,b_mu=1,
                        a_alpha=5,b_alpha=1,m_beta=0,s2_beta=1,
                        cs_zeta=5,cs_beta=.3,B=100000))
c(z$acc_zeta, z$acc_beta)

talpha <- tail(z$alpha,BB)
ttheta <- tail(z$theta,BB)
tbeta <- tail(z$beta,BB)
tzeta <- tail(z$zeta,BB)
tmu <- tail(z$mu,BB)

par(mfrow=c(4,1))
  plot(talpha,main=bquote(alpha),ylab='',type='l')
  plot(tzeta,main=bquote(zeta),ylab='',type='l')
  plot(tmu,main=bquote(mu),ylab='',type='l')
  plot(tbeta,main=bquote(beta),ylab='',type='l')
par(mfrow=c(1,1))

ut <-  apply(ttheta,1,function(x) length( unique(x) ))
tut <- table(ut)
plot(tut / sum(tut), main=paste("Mean Number of Clusters =",mean(ut)))

plot(ttheta[,1],cex=.1,col=rgb(.1,.1,.1,.2),main="theta")
for ( k in 1:ncol(z$theta) )
  points(ttheta[,k],cex=.1,col=rgb(.1,.1,.1,.2))

plot(tail(z$mu,BB), tail(z$zeta,BB),type='p',pch=20)

# Post Pred
y0.pred <- matrix(0,BB,n)
y.pred <- matrix(0,BB,n)
theta.pred <- apply(matrix(1:(BB)),1,function(b) {
            th <- 0

            if (talpha[b] / (talpha[b]+n) > runif(1)) {
              th <- rgamma(1,tzeta[b],rate=tzeta[b]/tmu[b])
            } else {
              th <- sample(ttheta[b,],1)
            }

            th
      })
for (i in 1:n) y0.pred[,i] <- rpois(BB, theta.pred*exp(cs(x)[i]*tbeta) )
for (i in 1:n) y.pred[,i] <- rpois(BB, ttheta[,i]*exp(cs(x)[i]*tbeta) )
ord <- order(x)
ymax <- max(apply(y0.pred,2,mean),apply(y.pred,2,mean),y)
plot(x[ord],y[ord],col="grey80",cex=2,pch=20,ylim=c(0,ymax))
lines(x[ord],apply(y0.pred,2,mean)[ord],col="green",pch=20,type='o')
points(x[ord],apply(y.pred,2,mean)[ord],col="blue",pch=20)

# Simple Poisson Regression ######
system.time( sourceCpp("regpois.cpp") )
system.time(z1 <- regpois(y,cs(x),zeta=.9368,mu=.39,m=0,s2=.1,cs_beta=.2,B=1000000))
z1$acc_beta
t1theta <- tail(z1$theta,BB)
t1beta <- tail(z1$beta,BB)
par(mfrow=c(2,1))
plot(t1theta,type='l',main=mean(t1theta))
plot(t1beta,type='l',main=mean(t1beta))
par(mfrow=c(1,1))

#plot(tail(rgp$beta,BB), tail(rgp$theta,BB), col=rgb(.5,.5,.5,.1:BB/BB),pch=20,cex=.1)
y1.pred <- apply(matrix(1:n), 1, 
                 function(i)
                 rpois( BB, t1theta*exp(t1beta*cs(x)[i]) ))
y1max <- max(apply(y1.pred,2,mean),apply(y1.pred,2,mean),y)
lines(x[ord],apply(y1.pred,2,mean)[ord],col="red",pch=20,type='o')


