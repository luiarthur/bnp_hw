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
viewYearJuly <- function(yr) {
  ind <- which(uyears==yr)
  quilt.plot(ylatlon[,2],ylatlon[,1],Y[ind,],
             fg='grey90',bty='n',main=yr,
             ylim=range(ylatlon[,1])+c(-1,1),
             xlim=range(ylatlon[,2])+c(-1,1),
             breaks=seq(14,40,len=101),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewYearJuly(1989) #1985 - 2004

#rr <- 1e4; cc <- 100
#XX <- matrix(sample(c(rnorm(3,0,1000)),rr*cc,repl=T),rr,cc)
#system.time(a <- unique(XX))
#system.time(b <- uniqueRows2(XX))
#all(a==b)
#uniqueRows(XX)

sourceCpp("sdp.cpp")

s_new <- matrix(0:50,ncol=2)
out <- sdp(Y, s_ner , D, beta_mu=0, beta_s2 = 100,
           tau2_a = 2, tau2_b = 10, alpha_a = 3, alpha_b=3,
           sig2_a = 2, sig2_b = 10, phi_a=0, phi_b=25, B=10)

