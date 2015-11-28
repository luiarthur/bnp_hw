library(Rcpp)
Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11, for RcppArmadillo
source("../../../R_Functions/plotPost.R",chdir=T)
cmaq <- read.csv("../../bnp_gp/data/CMAQ.csv")
o3 <- read.csv("../../bnp_gp/data/Ozone.csv")
predloc <- read.csv("../../bnp_gp/data/PredLocs.csv")
cs <- function(x) ( x - mean(x) ) / sd(x)
x <- dat$length
y <- dat$faults
n <- nrow(dat)


