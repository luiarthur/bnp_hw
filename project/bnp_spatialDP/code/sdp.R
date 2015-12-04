library(LatticeKrig) # quilt.plot
library(maps) # map()
library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
s_new <- read.csv("../data/predlocs.dat")
set.seed(1)
s_new1 <- s_new[sample(1:nrow(s_new),300,repl=T),]

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
ylatlon <- Yout[,1,3:4]
D <- as.matrix(dist(ylatlon))
# lon,lat,val
viewYearJuly <- function(yr,y) {
  ind <- which(uyears==yr)
  quilt.plot(ylatlon[,2],ylatlon[,1],y[ind,],
             fg='grey90',bty='n',main=yr,
             ylim=range(ylatlon[,1])+c(-1,1),
             xlim=range(ylatlon[,2])+c(-1,1),
             breaks=seq(14,40,len=101),
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

s_new <- matrix(1:50,ncol=2)
sourceCpp("sdp.cpp")
out <- sdp(Y, s_new , D, beta_mu=0, beta_s2 = 100,
           tau2_a = 2, tau2_b = 10, alpha_a = 1, alpha_b=1,
           sig2_a = 2, sig2_b = 5, phi_a=0, phi_b=.1, B=500)

par(mfrow=c(5,1),mar=c(0,4.5,1,2),fg='grey30',bty='l')
plot(out$beta[-1],type='l',xaxt='n')
par(mar=c(0,4.5,0,2))
plot(out$alpha[-1],type='l',xaxt='n')
plot(out$tau2[-1],type='l',xaxt='n') # about 11.25
plot(out$sig2[-1],type='l',xaxt='n') # Challenge
par(mar=c(3,4.5,0,2))
plot(out$phi[-1],type='l')  # Challenge
par(mfrow=c(1,1),mar=c(5,4,4,2)+.1)

ot <- out$theta
dim(ot)
B <- dim(ot)[3]
ot[,,B]
unique(ot[,,B]) # I should see clustering across times
viewYearJuly(1998,ot[,,B])
viewYearJuly(1998,Y)

mean( apply(Y,2,var) ) # empirical sig2. Difficult.
mean( apply(Y,1,var) ) # empirical tau


