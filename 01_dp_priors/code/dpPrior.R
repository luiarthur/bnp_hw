# Simulation of Dirichlet process prior realizations

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
dp <- function(N=1,a,pG,xlim=c(0,1),n=1000) {
  x <- seq(xlim[1],xlim[2],length=n)
  x <- sort(c(-1e1000,x))
  dG0 <- pG(x[-1]) - pG(x[-length(x)])
  out <- rdir(N,a*dG0)
  G <- t(apply(out,1,cumsum))
  list("G"=G, "x"= x[-1])
}

# DP using Sethuraman's construction
dp_stickbreak <- function(N=1,a,rG,xlim=c(0,1),J=1000) {
  x <- seq(xlim[1],xlim[2],length=J)

  z <- rbeta(J,1,a)
  th <- rG(J)
  w <- rep(0,J)
  w[1] <- z[1]
  for (l in 2:J) { # OPTIMIZE?
    w[l] <- z[l] * prod(1-z[1:(l-1)])
  }
  # v1
  #G <- sample(th,J,replace=T,prob=w)
  #G
 
  G <- matrix(0,N,J)
  for (i in 1:N) {
    for (l in 1:J) {
      G[i,l] <- sum(w[th <= x[l]])
    }
  }

  out <- list("G"=G, "x"=th)
  #plot(x,G,type="l")
  out
}
gx <- dp_stickbreak(a=1,rG=function(n) rnorm(n))

### Plotting Function
dp.post <- function(X,...) {
  plot(0,xlim=range(X$x),ylim=c(0,1),cex=0,...)
  for (i in 1:nrow(X$G)) {
    lines(X$x,X$G[i,],type="l",col=rgb(.4,.4,.4,.1))
  }
}
