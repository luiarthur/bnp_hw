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
             fg='grey90',bty='n',
             ylim=range(ylatlon[,1])+c(-1,1),
             xlim=range(ylatlon[,2])+c(-1,1),
             breaks=seq(14,40,len=101),
             col= colorRampPalette(c('dark blue','grey90','dark red'))(100))
  map('county',add=T,col='grey')
  map('state',add=T,col='grey60',lwd=2)
}
viewYearJuly(1989) #1985 - 2004

#XX <- matrix(sample(0:2,10*3,replace=T),10,3)
#unique(XX)
#uniqueRows(XX)

sourceCpp("sdp.cpp")

