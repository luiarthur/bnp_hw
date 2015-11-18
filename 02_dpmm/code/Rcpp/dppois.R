library(Rcpp)
dat <- read.table("../../dat/fabric.dat",header=T)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11
system.time( sourceCpp("dppois.cpp") )

z <- dppois(y,x,1,1,1,1,1,1,1,1,1,1,B=1000)
