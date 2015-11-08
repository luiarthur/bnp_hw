useRcpp <- F
y <- read.table("../../dat/hw2.dat")[[1]]

B <- 10000
burn <- B * .3

if (useRcpp) {
# cpp Gibbs:####################################
  library(Rcpp)
  sourceCpp("dp_arma.cpp")

  system.time(out <- gibbs(y,a=2,s=1,cs=3,B)) # 30 times faster than julia

  tab <- table(out[B,])
  mt <- apply(tail(out,B-burn),2,mean)
  mns <- apply(tail(out,B-burn),1,function(x) length(unique(x)))
  tabm <- table(mns) / sum(table(mns))
  plot(tabm)

  plot(y,pch=20,col='grey',cex=2)
  points(mt,pch=20,col=rgb(0,0,1,.5))
} else {
#R Gibbs: ######################################
  n <- length(y)
  theta <- matrix(0,B,n)
  for (b in 2:B) {
    for (i in 1:n) {
      tb <- theta[b,]

    }
  }
}
