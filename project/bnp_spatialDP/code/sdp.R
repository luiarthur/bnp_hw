library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
s_new <- read.csv("../data/predlocs.dat")
set.seed(1)

load("../data/Y.RData")

Ydim <- dim(Yout)
locs <- Ydim[1]
datcols <- Ydim[3]
ntimes <- Ydim[2]

uyears <- unique(Yout[1,,1])
TT <- length(uyears)
july <- which((1:ntimes)%%12%%7==0 & (1:ntimes)%%12!=0)
Y <- t(Yout[,july,5]) # Y is 20 x 100
                      # 20 years of july max temperatures
                      # 100 stations
#Y <- scale(Y,center=TRUE,scale=FALSE)[1:nrow(Y),]
ylatlon <- Yout[,1,3:4]
D <- as.matrix(dist(ylatlon))
# lon,lat,val
viewPred <- function(x,latlon,main.plot='',bks=c(14,40)) {
  quilt.plot(latlon[,2],latlon[,1],x,
             fg='grey90',bty='n',main=main.plot,
             ylim=range(latlon[,1])+c(-1,1),
             xlim=range(latlon[,2])+c(-1,1),
             breaks=seq(bks[1],bks[2],len=length(x)+1),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(length(x)))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}


viewYearJuly <- function(yr,y,bks=seq(14,40,len=101)) {
  ind <- which(uyears==yr)
  quilt.plot(ylatlon[,2],ylatlon[,1],y[ind,],
             fg='grey90',bty='n',main=yr,
             ylim=range(ylatlon[,1])+c(-1,1),
             xlim=range(ylatlon[,2])+c(-1,1),
             breaks=bks,
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewYearJuly(1989,Y) #1985 - 2004


#rr <- 100000; cc <- 3
#XX <- matrix(sample(c(-100,0,100),rr*cc,repl=T),rr,cc)
#system.time(a <- unique(XX))
#system.time(e <- uniqueRows(XX))
#all(a == e)

sourceCpp("sdp.cpp")
system.time(
#out <- sdp(t(Yout[,,5]), s_new , D, beta_mu=0, beta_s2 = 1,
out <- sdp(Y, D, beta_m=30, beta_s2 = 100,
           tau2_a = 2, tau2_b = 100, alpha_a = 1, alpha_b=.0001,
           sig2_a = 2, sig2_b = .5, phi_a=.3, phi_b=.3, L=10, B=50000)
)

burn <- round(B*.5)
par(mfrow=c(6,1),mar=c(0,4.5,1,2),fg='grey30',bty='l')
plot(out$beta[-c(1:burn)],type='l',xaxt='n',ylab=bquote(beta))
plot(out$alpha[-c(1:burn)],type='l',xaxt='n',ylab=bquote(alpha))
par(mar=c(0,4.5,0,2))
plot(out$tau2[-c(1:burn)],type='l',xaxt='n',ylab=bquote(tau^2)) # about 11.25
plot(out$sig2[-c(1:burn)],type='l',xaxt='n',ylab=bquote(sigma^2)) # about 80
plot(out$theta[20,62,-c(1:burn)],type='l',xaxt='n',ylab=bquote(theta))
par(mar=c(3,4.5,0,2))
plot(out$phi[-c(1:burn)],type='l',ylab=bquote(phi)) #about .05
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1)

ot <- out$theta
dim(ot)
#ot[,,7]
B <- dim(ot)[3]
b <- B
uB <- uniqueRows(ot[,,b]) # I should see clustering across times
unique(uB); nrow(uB)
nURows <- apply(ot[,,round(B*.8):B],3,function(x) nrow(uniqueRows(x)))
plot(table(nURows)/sum(table(nURows)))
apply(uB,1,function(x) length(matchRows(ot[,,b], x)))
ind <- matchRows(ot[,,b], uB[1,]); ind # Change the index to see what happens
par(mfrow=c(ifelse(length(ind)>1,ceiling(length(ind)/2),1),ifelse(length(ind)>1,2,1)),
    mar=c(1,4.5,1,2),fg='grey30',bty='l')
for (ii in ind) {
  viewYearJuly(1985+ii,Y)
  #viewPred(t(Yout[,ii+1,5]),ylatlon,bks=range(Yout[,ind+1,5]))
}
par(mfrow=c(1,1))

mt <- t(apply(ot,1,function(x) apply(x,1,mean))) # mean theta 20 x 100
par(mfrow=c(1,3))
viewYearJuly(1989,Y)
viewYearJuly(1989,mt)
viewYearJuly(1989,mt-Y,bks=seq(-2,2,len=dim(ot)[2]+1))
par(mfrow=c(1,1))

mean( apply(Y,2,var) ) # empirical sig2. Difficult.
mean( apply(Y,1,var) ) # empirical tau

### Predictions at Observed Locations:
# Obtain theta0
keep <- round(B*.3)
post.a <- tail(out$alpha,keep)
post.b <- tail(out$beta,keep)
post.p <- tail(out$phi,keep)
post.s2 <- tail(out$sig2,keep)
post.t2 <- tail(out$tau2,keep)
post.th <- out$theta[,,(B-keep+1):B]
pred.theta0 <- matrix(0,keep,dim(ot)[2])
for (bb in 1:keep) {
  idx_new <- sample(0:(dim(ot)[1]), 1, prob=c(post.a[bb], rep(1,dim(ot)[1]) ))
  if (idx_new == 0) {
    pred.theta0[bb,] <- mvrnorm(rep(0,dim(ot)[2]),post.s2[bb]*Hn(post.p[bb],D))
  } else {
    pred.theta0[bb,] <- post.th[idx_new,,bb]
  }
  cat("\r",bb)
}
round(apply(pred.theta0,2,mean),4)

pred.y0 <- t(apply(matrix(1:keep),1,function(i) 
                   mvrnorm(pred.theta0[i,]+post.b[i],
                           post.t2[i]*diag(ncol(Y)))))
apply(Y,2,mean)
apply(pred.y0,2,mean)
apply(pred.y0,2,var)
par(mfrow=c(1,2))
viewPred(apply(pred.y0,2,var),ylatlon, "Posterior Predictive Variance",bks=c(0,40))
viewPred(apply(Y,2,var), ylatlon, "Variance of Data",bks=c(0,3))
par(mfrow=c(1,1))

#hist(Y[,3],prob=TRUE)
#lines(density(pred.y0[,3]),col='red')


# Print this:
par(mfrow=c(1,2))
  viewPred(apply(pred.y0,2,mean), ylatlon, "Posterior Predictive Median")
  viewPred(apply(Y,2,mean), ylatlon, "Data Median")
par(mfrow=c(1,1))

### Predictions at NEW Locations:

nn <- ncol(Y)
n_new <- 1000
s0 <- s_new[sample(1:nrow(s_new),n_new,replace=F),]
D_new <- as.matrix(dist(s0))
pred.new.theta <- matrix(0,keep,n_new)
pred.new.y <- matrix(0,keep,n_new)

for (bb in 1:keep) {
  pred.new.theta[bb,] <- mvrnorm(rep(0,n_new),
                                post.s2[bb]*Hn(post.p[bb],D_new))
  pred.new.y[bb,] <- mvrnorm(pred.new.theta[bb,]+post.b[bb],
                             post.t2[bb]*diag(n_new));
  cat("\r",bb)
}

apply(pred.new.theta,2,mean)
apply(pred.new.theta,2,var)
apply(pred.new.y,2,mean)
apply(pred.new.y,2,var)

par(mfrow=c(1,3))
  viewPred(apply(pred.y0,2,mean), ylatlon, "Posterior Predictive Median")
  viewPred(apply(Y,2,mean), ylatlon, "Data Median")
  viewPred(apply(pred.new.y,2,mean)[1:n_new], s0[,2:1], "Posterior Predictive Mean",bks=c(0,40))
par(mfrow=c(1,1))

