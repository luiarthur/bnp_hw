library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
load("../data//Y.RData")
Y <- Yout; rm(Yout)

sourceCpp("sdp.cpp")

