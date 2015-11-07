library(Rcpp)
sourceCpp("dp_arma.cpp")
system.time(x <- rdir(100000,1:100))

y1 <- read.table("../../dat/hw2.dat")[[1]]
plot(density(y1))


