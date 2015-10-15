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

dp_stickbreak <- function(N=1,a,rG,xlim=c(0,1), J=NULL, eps=.01, printProg=T, K=100) {
  if (is.null(J)) {
    J <- log(eps) / (log(a) - log(a+1))
    J <- round(J)
    print(J)
  }

  #x <- seq(xlim[1],xlim[2],length=J) # Change J to something else.
  x <- seq(xlim[1],xlim[2],length=K)
  Z <- matrix(rbeta(J*N,1,a),N,J)
  th <- matrix(rG(J*N),N,J)
  W <- matrix(0,N,J)
  #G <- matrix(0,N,J)
  G <- matrix(0,N,K)
  W[,1] <- Z[,1]

  for (i in 1:N) {
    for (l in 1:J) { # OPTIMIZE?
      if (l>1) W[i,l] <- Z[i,l] * prod(1-Z[i,1:(l-1)])
    }
    for (k in 1:K) { # OPTIMIZE?
      G[i,k] <- sum(W[i,th[i,] <= x[k]])
      it <- (i-1)*K+l
      if (printProg) cat("\r Progress: ", it/(N*K) )
    }
  }

  out <- list("G"=G, "x"=x, "J"=J)
  out
}

mdp <- function(N=1,rA,pG,xlim=c(0,1),n=100) {
  x <- seq(xlim[1],xlim[2],length=n)
  x <- c(-1e1000,x)
  a <- rA(N) # prior for alpha
  adG0 <- apply(matrix(a),1, function(aa) aa*pG(x[-1]) - aa*pG(x[-length(x)]))

  y <- matrix(rgamma(n*N, adG0, 1), N, n, byrow=T)
  rowsums <- y %*% rep(1,n)
  out <- y / as.vector(rowsums)

  G <- t(apply(out,1,cumsum))
  list("G"=G, "x"= x[-1],"a"=a)
}

dp.post <- function(X,col.lines=rgb(.4,.4,.4,.1),xlim.def=range(X$x),...) {
  plot(0,xlim=xlim.def,ylim=c(0,1),cex=0,bty="n",las=1,
       col.axis=rgb(.3,.3,.3),fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),
       col.main=rgb(.3,.3,.5),...)
  xx <- X$x
  N <- nrow(X$G)
  for (i in 1:N) {
    if (i > N - 3) {
      #lines(xx,X$G[i,],type="l",col=rgb(.2,.2,.2))
      lines(xx,X$G[i,],type="l",col="black")
    } else {
      lines(xx,X$G[i,],type="l",col=col.lines)
    }
  }
}


# dp(N=1000, a=1, pG = function(x) pnorm(x), xlim=c(-3,3), n=100)
# dp_stickbreak(N=1000, a=1, rG = function(n) rnorm(n), xlim=c(-3,3), K=100)
# mdp(N=1000, rA=function(n) rgamma(n,1,2), pG=function(x) pnorm(x), xlim=c(-3,3))
