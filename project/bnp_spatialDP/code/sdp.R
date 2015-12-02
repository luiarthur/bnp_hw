library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
load("../data//Y.RData")

Ydim <- dim(Yout)
locs <- Ydim[1]
datcols <- Ydim[3]
ntimes <- Ydim[2]

july <- which((1:ntimes)%%12%%7==0 & (1:ntimes)%%12!=0)
Y <- t(Yout[,july,5]) # Y is 20 x 100
                      # 20 years of july max temperatures
                      # 100 stations
ylatlon <- Yout[,1,3:4]


#XX <- matrix(sample(0:2,10*3,replace=T),10,3)
#unique(XX)
#uniqueRows(XX)

sourceCpp("sdp.cpp")
