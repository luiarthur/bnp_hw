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
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo

BB <- 5000
# DP ###
system.time( sourceCpp("dppois.cpp") )#, showOutput=T) )
system.time(m3 <- dppois(y,cs(x),a_zeta=1,b_zeta=1,a_mu=2,b_mu=1,
                        a_alpha=1,b_alpha=1,m_beta=0,s2_beta=1,
                        cs_zeta=3,cs_beta=.2,B=1000000))
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
  plot(x[ord],y[ord],col="grey80",cex=2,pch=20,ylim=c(0,m3_ymax),bty='l',fg='grey50',main="Posterior Predictive",ylab="y",xlab="x")
  #lines(x[ord],apply(m3_y0.pred,2,mean)[ord],col="green",pch=20,type='o')
  points(x[ord],apply(m3_y.pred,2,mean)[ord],col="gold",pch=20,cex=1.5)
}
m3.post.pred()


# Simple Poisson Regression ######
system.time( sourceCpp("regpois.cpp") )
system.time(m1 <- regpois(y,cs(x),zeta=1,mu=1,m=0,s2=1,cs_beta=.2,B=1000000))
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
lines(x[ord],apply(m1_y0.pred,2,mean)[ord],col="red",pch=20,type='o',lwd=2)

plot(density(apply(m3_theta,1,mean)),col='red'); lines(density(m1_theta))
plot(density(m3_beta),col='red'); lines(density(m1_beta))

# Hierarchical Model
system.time( sourceCpp("hier.cpp") )
system.time(m2 <- hier(y,cs(x),a_zeta=1,b_zeta=1,a_mu=2,b_mu=1,
                       m_beta=0,s2_beta=1,cs_zeta=3,cs_beta=.22,B=1000000))
c(m2$acc_zeta, m2$acc_beta)

m2_theta <- tail(m2$theta,BB)
m2_beta <- tail(m2$beta,BB)
m2_zeta <- tail(m2$zeta,BB)
m2_mu <- tail(m2$mu,BB)

par(mfrow=c(2,2))
  plot.post(m2_zeta,main=bquote(zeta),ylab='')
  plot.post(m2_mu,main=bquote(mu),ylab='')
  plot.post(m2_beta,main=bquote(beta),ylab='')
  plot(m2_theta[,1],cex=.1,col=rgb(.1,.1,.1,.2),main=bquote(theta),fg='grey30',bty='l')
  for ( k in 1:ncol(m2_theta) )
    points(m2_theta[,k],cex=.1,col=rgb(.1,.1,.1,.2))
par(mfrow=c(1,1))
dev.off()

m2_y.pred <- matrix(0,BB,n)
for (i in 1:n) m2_y.pred[,i] <- rpois(BB, m2_theta[,i]*exp(cs(x)[i]*m2_beta) )



####
pdf("../../latex/img/poisPostPred.pdf")
m3.post.pred()
points(x[ord],apply(m2_y.pred,2,mean)[ord],col="blue",pch=20,cex=1.5)
lines(x[ord],apply(m1_y0.pred,2,mean)[ord],col="red",pch=20,type='o',lwd=2)
legend("topleft",legend=c("Simple Poisson Regression","Hierarchical Poisson Regression","DP Mixture Regression","Data"),bty='n',
       col=c("red","blue","gold","grey"),lwd=3)
dev.off()

plot(density(apply(m3_theta,1,mean)),col='red'); lines(density(m1_theta)); lines(density(apply(m2_theta,1,mean)),col='blue');
plot(density(m3_beta),col='red'); lines(density(m1_beta)); lines(density(m2_beta),col='blue');

#CPO: ##############
cpo <- function(post_theta,post_beta) {
  M <- matrix(0,BB,n)
  for (i in 1:n) {
    M[,i] <- 1 / dpois(y[i],post_theta[,i] * exp(cs(x)[i]*post_beta) )
  }
  denom <- apply(M,2,mean)
  1 / denom
}
(m1_cpo <- cpo(matrix(rep(m1_theta,n),BB), m1_beta))
(m2_cpo <- cpo(m2_theta, m2_beta))
(m3_cpo <- cpo(m3_theta, m3_beta))

prod(m3_cpo / m2_cpo)
prod(m3_cpo / m1_cpo)
prod(m2_cpo / m1_cpo)

pdf("../../latex/img/cpo.pdf")
   plot(x[ord], log(m2_cpo / m1_cpo)[ord], type='p', pch=20, cex=1.5, col='red', bty='l',fg='grey',xlab='x',ylab='Log cpo ratio')
  lines(x[ord], log(m3_cpo / m1_cpo)[ord], type='p', pch=20, cex=1.5, col='green')
  lines(x[ord], log(m3_cpo / m2_cpo)[ord], type='p', pch=20, cex=1.5, col='blue')
  points(x[ord],log(y[ord]),col='grey',pch=20,cex=1)
  legend("topleft",legend=c("log(cpo3 / cpo2)", "log(cpo3 / cpo1)","log(cpo2 / cpo1)","log(y)"),col=c("blue","green","red","grey"),lwd=2,bty='n')
dev.off()

