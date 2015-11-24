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

BB <- 5000
# DP ###
system.time( sourceCpp("dppois.cpp") )#, showOutput=T) )
system.time(m3 <- dppois(y,cs(x),a_zeta=1,b_zeta=1,a_mu=1,b_mu=1,
                        a_alpha=1e-10,b_alpha=1,m_beta=0,s2_beta=1,
                        cs_zeta=9,cs_beta=.2,B=100000))
c(m3$acc_zeta, m3$acc_beta)

m3_alpha <- tail(m3$alpha,BB)
m3_theta <- tail(m3$theta,BB)
m3_beta <- tail(m3$beta,BB)
m3_zeta <- tail(m3$zeta,BB)
m3_mu <- tail(m3$mu,BB)

par(mfrow=c(2,2))
  plot.post(m3_alpha,main=bquote(alpha),ylab='')
  plot.post(m3_zeta,main=bquote(zeta),ylab='')
  plot.post(m3_mu,main=bquote(mu),ylab='')
  plot.post(m3_beta,main=bquote(beta),ylab='')
par(mfrow=c(1,1))
dev.off()

plot(m3_mu, m3_zeta,type='p',pch=20)

ut <-  apply(m3_theta,1,function(x) length( unique(x) ))
tut <- table(ut)
plot(tut / sum(tut), main=paste("Mean Number of Clusters =",mean(ut)))

plot(m3_theta[,1],cex=.1,col=rgb(.1,.1,.1,.2),main="theta")
for ( k in 1:ncol(m3_theta) )
  points(m3_theta[,k],cex=.1,col=rgb(.1,.1,.1,.2))


# Post Pred
m3_y0.pred <- matrix(0,BB,n)
m3_y.pred <- matrix(0,BB,n)
m3_theta.pred <- apply(matrix(1:(BB)),1,function(b) {
            th <- 0

            if (m3_alpha[b] / (m3_alpha[b]+n) > runif(1)) {
              th <- rgamma(1,m3_zeta[b],rate=m3_zeta[b]/m3_mu[b])
            } else {
              th <- sample(m3_theta[b,],1)
            }

            th
      })
for (i in 1:n) m3_y0.pred[,i] <- rpois(BB, m3_theta.pred*exp(cs(x)[i]*m3_beta) )
for (i in 1:n) m3_y.pred[,i] <- rpois(BB, m3_theta[,i]*exp(cs(x)[i]*m3_beta) )
ord <- order(x)
m3_ymax <- max(apply(m3_y0.pred,2,mean),apply(m3_y.pred,2,mean),y)

m3.post.pred <- function() {
‘"l"’, ‘"7"’, ‘"c"’, ‘"u"’, or ‘"]"’
  plot(x[ord],y[ord],col="grey80",cex=2,pch=20,ylim=c(0,m3_ymax),bty='l',fg='grey50')
  lines(x[ord],apply(m3_y0.pred,2,mean)[ord],col="green",pch=20,type='o')
  points(x[ord],apply(m3_y.pred,2,mean)[ord],col="darkgreen",pch=20)
}
m3.post.pred()


# Simple Poisson Regression ######
system.time( sourceCpp("regpois.cpp") )
system.time(m1 <- regpois(y,cs(x),zeta=6.6,mu=8.29,m=0,s2=1,cs_beta=.2,B=1000000))
m1$acc_beta
m1_theta <- tail(m1$theta,BB)
m1_beta <- tail(m1$beta,BB)
par(mfrow=c(2,1))
plot(m1_theta,type='l',main=mean(m1_theta))
plot(m1_beta,type='l',main=mean(m1_beta))
par(mfrow=c(1,1))

#plot(tail(rgp$beta,BB), tail(rgp$theta,BB), col=rgb(.5,.5,.5,.1:BB/BB),pch=20,cex=.1)
m1_y0.pred <- apply(matrix(1:n), 1, function(i) rpois( BB, m1_theta*exp(m1_beta*cs(x)[i]) ))

m3.post.pred()
lines(x[ord],apply(m1_y0.pred,2,mean)[ord],col="red",pch=20,type='o')

plot(density(apply(m3_theta,1,mean)),col='red'); lines(density(m1_theta))
plot(density(m3_beta),col='red'); lines(density(m1_beta))


