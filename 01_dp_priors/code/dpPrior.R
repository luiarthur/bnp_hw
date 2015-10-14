# Simulation of Dirichlet process prior realizations
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")

color.btwn <- function(x,ylo,yhi,from,to,col.area="grey") {
  x <- c(x,rev(x))
  y <- c(yhi,rev(ylo))

  polygon(c(x[x>=from & x<= to]),
          c(y[x>=from & x<=to]),
          col=col.area,border=F)
}

rdir <- function(N,a) {
  k <- length(a)
  x <- matrix(rgamma(k*N, a, 1), N, k, byrow=T)
  rowsums <- x %*% rep(1,k)
  x / as.vector(rowsums)
}

# Using the dirichlet distribution (Ferguson)
dp <- function(N=1,a,pG,xlim=c(0,1),n=100) {
  x <- seq(xlim[1],xlim[2],length=n)
  x <- sort(c(-1e1000,x))
  dG0 <- pG(x[-1]) - pG(x[-length(x)])
  out <- rdir(N,a*dG0)
  G <- t(apply(out,1,cumsum))
  list("G"=G, "x"= x[-1])
}

### Plotting Function
dp.post <- function(X,col.lines=rgb(.4,.4,.4,.1),xlim.def=range(X$x),...) {
  plot(0,xlim=xlim.def,ylim=c(0,1),cex=0,bty="n",las=1,
       col.axis=rgb(.3,.3,.3),fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),
       col.main=rgb(.3,.3,.5),...)
  xx <- X$x
  print(nrow(X$G))
  for (i in 1:nrow(X$G)) {
    if (i<=3) {
      lines(xx,X$G[i,],type="l",col=rgb(.2,.2,.2))
    } else {
      lines(xx,X$G[i,],type="l",col=col.lines)
    }
  }

  lines(xx,apply(X$G,2,mean),col="blue",lwd=2)
  lines(xx,apply(X$G,2,var),col="red",lwd=2)
  minor <- function() 
    plot(xx,apply(X$G,2,var),col="red",lwd=1,bty="n",type="l",
         col.axis=rgb(.3,.3,.3),fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),
         col.main=rgb(.3,.3,.5))
  plot.in.plot(minor,"topleft")
}

# DP using Sethuraman's construction
dp_stickbreak <- function(N=1,a,rG,xlim=c(0,1), J=NULL, eps=.01, printProg=T) {
  if (is.null(J)) {
    J <- log(eps) / (log(a) - log(a+1))
    J <- round(J)
    print(J)
  }

  x <- seq(xlim[1],xlim[2],length=J)
  Z <- matrix(rbeta(J*N,1,a),N,J)
  th <- matrix(rG(J*N),N,J)
  W <- matrix(0,N,J)
  G <- matrix(0,N,J)
  W[,1] <- Z[,1]

  for (i in 1:N) {
    for (l in 1:J) { # OPTIMIZE?
      if (l>1) W[i,l] <- Z[i,l] * prod(1-Z[i,1:(l-1)])
      G[i,l] <- sum(W[i,th[i,] <= x[l]])
      it <- (i-1)*J+l
      if (printProg) cat("\r Progress: ", it/(N*J) )
    }
  }

  out <- list("G"=G, "x"=x, "J"=J)
  out
}
# Example
# gx <- dp_stickbreak(N=1000, a=3,rG=function(n) rnorm(n), xlim=c(-3,3), eps=1e-4)
# dp.post(gx,col.lines=rgb(.4,.4,.4,.05),ylab="F(x)",xlab="x",main="DP")


# 2a
# Setharuman
par(mfrow=c(1,3))
avec <- c(1,10,100)
for (a in avec) {
  gs <- dp_stickbreak(N=1000, a=a,rG=function(n) rnorm(n), xlim=c(-3,3), eps=1e-4)
  #main <- bquote(paste("DP(",.(al),"G"[0],")"))
  dp.post(gs,col.lines=rgb(.4,.4,.4,.05),ylab="F(x)",xlab="x",
          main=bquote("G ~ DP("~.(a)~","~G[0]~")"~" - Sethuraman's Construction"))
}
par(mfrow=c(1,1))

# Ferguson
par(mfrow=c(1,3))
avec <- c(1,10,100)
for (a in avec) {
  gs <- dp(N=1000, a=a,pG=function(n) pnorm(n), xlim=c(-3,3))
  #main <- bquote(paste("DP(",.(al),"G"[0],")"))
  dp.post(gs,col.lines=rgb(.4,.4,.4,.05),ylab="F(x)",xlab="x",
          main=bquote("G ~ DP("~.(a)~","~G[0]~")"~" - Ferguson's Construction"))
}
par(mfrow=c(1,1))

#source("dpPrior.R")
