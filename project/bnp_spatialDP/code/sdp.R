library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
load("../data//Y.RData")

#XX <- matrix(sample(0:2,10*3,replace=T),10,3)
#unique(XX)
#uniqueRows(XX)


sourceCpp("sdp.cpp")

