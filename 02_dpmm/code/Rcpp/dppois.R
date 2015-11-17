library(Rcpp)
dat <- read.table("../../dat/fabric.dat",header=T)
x <- dat$length
y <- dat$faults
n <- nrow(dat)
#plot(x,y)

sourceCpp("dp_arma.cpp")
system.time(x <- rdir(100000,1:100))


