library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
library(doMC)
registerDoMC(as.numeric(system("nproc",intern=TRUE))/2)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=TRUE)
s_new <- read.csv("../data/predlocs.dat")
set.seed(1)

load("../data/Y.RData")

Ydim <- dim(Yout)
s <- Yout[,1,3:4] # latlon
colnames(s) <- c("lat","lon")
n <- nrow(s)
uyears <- unique(Yout[1,,1])
T <- length(uyears)
july <- which((1:(Ydim[2]))%%12%%7==0 & (1:(Ydim[2]))%%12!=0)
Y <- t(Yout[,july,5]) # Y is 20 x 100
                      # 20 years of july max temperatures
                      # 100 stations
D <- as.matrix(dist(s))
viewPred <- function(x,latlon,main.plot='',bks=range(x)) {
  quilt.plot(latlon[,2],latlon[,1],x,
             fg='grey90',bty='n',main=main.plot,
             ylim=range(latlon[,1])+c(-1,1),
             xlim=range(latlon[,2])+c(-1,1),
             breaks=seq(bks[1],bks[2],len=length(x)+1),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(length(x)))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
# Print this:
viewPred(Y[1,],s,bks=c(15,40))

mean( apply(Y,2,var) ) # empirical sig2. Difficult.
mean( apply(Y,1,var) ) # empirical tau

plotPosts <- function(OUT) {
  OT <- OUT$theta
  B <- dim(OT)[3]
  burn <- round(B*.8)
  par(mfrow=c(7,1),mar=c(0,4.5,1,2),fg='grey30',bty='l')
  plot(OUT$beta[-c(1:burn)],type='l',xaxt='n',ylab=bquote(beta))
  plot(OUT$alpha[-c(1:burn)],type='l',xaxt='n',ylab=bquote(alpha))
  par(mar=c(0,4.5,0,2))
  plot(OUT$tau2[-c(1:burn)],type='l',xaxt='n',ylab=bquote(tau^2))# ~13
  plot(OUT$sig2[-c(1:burn)],type='l',xaxt='n',ylab=bquote(sigma^2))#~100
  plot(OUT$the[20,39,-c(1:burn)],type='l',xaxt='n',yla=bquote(theta[39]))
  plot(OUT$the[20,62,-c(1:burn)],type='l',xaxt='n',yla=bquote(theta[62]))
  par(mar=c(3,4.5,0,2))
  plot(OUT$phi[-c(1:burn)],type='l',ylab=bquote(phi)) #about .05
  par(mfrow=c(1,1),mar=c(5,4,4,2)+.1)
}

y.x <- function(x,y_idx,x_idx,m,S) {
  mu.y <- m[y_idx]
  mu.x <- m[x_idx]
  S.yx <- S[y_idx,x_idx]
  S.xy <- t(S.yx)
  S.yy <- S[y_idx,y_idx]
  Sxx.i <- solve(S[x_idx,x_idx])
  E.y.x <- mu.y + S.yx %*% Sxx.i %*% (x-mu.x)
  #V.y.x <- S.yy - S.yx %*% Sxx.i %*% S.xy
  #list("m"=E.y.x, "S"=V.y.x)
  list("m"=E.y.x)
}
#########################################################################
sourceCpp("sdp.cpp")
system.time(
out <- sdp(Y, D, beta_m=30, beta_s2 = 100,
           tau2_a = 2, tau2_b = 13, alpha_a = 1, alpha_b=10,
           sig2_a = 2, sig2_b = 100, phi_a=1, phi_b=1, L=10, B=2500)
)
#update_tau2(2,100,100,20,Y,matrix(rnorm(2000,0,1),20),29)
ot <- out$theta
B <- length(out$alpha)
plotPosts(out)

### Predictions at Observed Locations:
keep <- round(B*.1)
post.a <- tail(out$alpha,keep)
post.b <- tail(out$beta,keep)
post.p <- tail(out$phi,keep)
post.s2 <- tail(out$sig2,keep)
post.t2 <- tail(out$tau2,keep)
post.th <- out$theta[,,(B-keep+1):B]

### New Locations:
m <- nrow(s_new)#2000 # number of new locations
s0 <- as.matrix(s_new[sample(1:nrow(s_new),m,replace=F),])[,2:1] 
colnames(s0) <- c("lat","lon")
D.mn <- as.matrix(dist(rbind(s0,s)))

post.th0 <- array(0,dim=c(T,m,keep))
dim(post.th0)
for (b in 1:keep) {
  th <- post.th[,,b]
  ur <- uniqueRows(th)
  J <- nrow(ur)
  for (j in 1:J) {
    idx <- which( apply(matrix(1:T),1,function(t) all(th[t,] == ur[j,]) ))
    lmS <- y.x(c(ur),1:m,(m+1):(m+n),rep(0,m+n),Hn(post.p[b],D.mn))
    for (i in idx) post.th0[i,,b] <- lmS$m
  }
  cat("\r",b,"/",keep)
}

pred.t0 <- matrix(0,keep,m+n)
for (b in 1:keep) {
  idx_new <- sample(0:T, 1, prob=c(post.a[b], rep(1,T) ))
  if (idx_new == 0) {
    pred.t0[b,] <- mvrnorm(rep(0,m+n),post.s2[b]*Hn(post.p[b],D.mn))
  } else {
    pred.t0[b,(m+1):(m+n)] <- post.th[idx_new,,b]
    pred.t0[b,1:m] <- post.th0[idx_new,,b]
  }
  cat("\r",b)
}
round(apply(pred.t0,2,mean),4)



get.pred.y0 <- function(i) {
  cat("\r",i/keep)
  c(mvrnorm(pred.t0[i,]+post.b[i], post.t2[i]*diag(m+n)))
}
# Sequential:
# pred.y0 <- t(apply(matrix(1:keep),1, get.pred.y0(i)))
# Parallel:
pred.y0 <- foreach(i=1:keep,.combine=rbind) %dopar% get.pred.y0(i)
apply(Y,2,mean)
apply(pred.y0,2,mean)
apply(pred.y0,2,var)
# Print this:
par(mfrow=c(1,2))
  viewPred(apply(Y,2,mean), s, "Data Mean")
  viewPred(apply(pred.y0,2,mean), rbind(s0,s), "Posterior Predictive Mean")
par(mfrow=c(1,1))
par(mfrow=c(1,2))
  viewPred(apply(Y,2,var), s, "Variance of Data",bks=c(0,20))
  viewPred(apply(pred.y0,2,var),rbind(s0,s), "Posterior Predictive Variance",bks=c(0,20))
par(mfrow=c(1,1))
