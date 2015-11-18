library(Rcpp)
dat <- read.table("../../dat/fabric.dat",header=T)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)

Sys.setenv("PKG_CXXFLAGS"="-std=c++11") # enable c++11
system.time( sourceCpp("dppois.cpp", showOutput=T) )
system.time(z <- dppois(y,x,1,1,1,1,1,1,1,1,100,.0001,B=10000))

plot(tail(z$alpha,2000),type='l')
plot(tail(z$mu,2000),type='l')
plot(tail(z$tau,2000),type='l')
plot(tail(z$beta,2000),type='l')


