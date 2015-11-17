library(Rcpp)
dat <- read.table("../../dat/fabric.dat",header=T)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)

system.time( sourceCpp("dppois.cpp") )

