rdir <- function(N,a) {
  k <- length(a)
  x <- matrix(rgamma(k*N, a, 1), N, k, byrow=T)
  rowsums <- x %*% rep(1,k)
  x / as.vector(rowsums)
}

# Using the dirichlet distribution (Ferguson)
dp <- function(N=1,a,pG,xlim=c(0,1),J=100) {
  if (length(xlim) <= 2) {
    x <- seq(xlim[1],xlim[2],length=J)
  } else {
    x <- xlim
  }

  dG0 <- c(pG(x[1]), pG(x[-1]) - pG(x[-length(x)]))
  out <- rdir(N,a*dG0)
  G <- t(apply(out,1,cumsum))
  list("G"=G, "x"= x)
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
      lines(xx,X$G[i,],type="s",col="black")
    } else {
      lines(xx,X$G[i,],type="s",col=col.lines)
    }
  }
}

# dp(N=1000, a=1, pG = function(x) pnorm(x), xlim=c(-3,3), n=100)
# dp_stickbreak(N=1000, a=1, rG = function(n) rnorm(n), xlim=c(-3,3), K=100)
# mdp(N=1000, rA=function(n) rgamma(n,1,2), pG=function(x) pnorm(x), xlim=c(-3,3))

dp.post.ci <- function(G, col.ci=rgb(.2,.2,.2,.5),ylim.def=c(0,1),xlim.def=range(G$x),lwd.EG=2,type.EG="l",
                       pch.EG=1,cex.EG=1,EG.col='blue',density.est=F,lwd.ci=1,...) {

  plot(0,cex=0,ylim=ylim.def,xlim=xlim.def,
       bty="n",las=1, col.axis=rgb(.3,.3,.3),
       fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),col.main=rgb(.3,.3,.4),...)
  
  g <- NULL
  if (!(density.est)) { # Plot cdf
    g <- G$G
  } else { # Plot pmf
    g <- t(apply(G$G,1,function(x) c(x[1], x[-1]-x[-length(x)]) ))
  }

  EG <- apply(g,2,function(x) mean(x,na.rm=T))
  VG <- apply(g,2,function(x) var(x,na.rm=T))
  qG <- apply(g,2,function(x) quantile(x,c(.025,.975),na.rm=T))
  glo <- qG[1,]
  ghi <- qG[2,]

  if (!(density.est)) {# cdf 95% CI
    lines(G$x,EG,col=EG.col,lwd=lwd.EG,type=type.EG,pch=pch.EG,cex=cex.EG) # E[G|y]
    color.btwn(G$x,glo,ghi,-100,100,col.area=col.ci)
  } else { # pmf 95% CI
    lines(G$x,EG,col=EG.col,lwd=lwd.EG,type=type.EG,pch=pch.EG,cex=cex.EG) # E[G|y]
    segments(G$x,glo,G$x,ghi,col=col.ci,lwd=lwd.ci)
  }
}

plot.cdf <- function(x,add=F,printProgress=F,...) {
  ux <- sort(unique(x))
  lux <- length(ux)
  cdf <- matrix(0,lux,2)

  for (i in 1:lux) {
    cdf[i,1] <- ux[i]
    if (i>1) cdf[i-1,2] <- sum(ux[i] > x)
    if (printProgress) cat("\r Progress: ",i,"/",lux)
  }
  if (printProgress) cat("\n")

  cdf[,2] <- cdf[,2] / length(x)
  cdf[lux,2] <- 1

  if (add) {
    lines(cdf[,1],cdf[,2],...)
  } else {
    plot(cdf[,1],cdf[,2],...)
  }
}

plot.pmf <- function(x,add=F,...) {
  x.table <- table(x)
  x.uniq <- as.numeric(names(x.table))
  x.counts <- as.numeric(x.table)

  lux <- length(x.uniq)
  pmf <- matrix(0,lux,2)

  x.prop <- x.counts / sum(x.counts)
  if (add) {
    lines(x.uniq, x.prop,...)
  } else {
    plot(x.uniq, x.prop,...)
  }

  out <- cbind(x.uniq, x.prop)
  colnames(out) <- c("x","prop")

  out
}
# Example:
#plot.cdf(rnorm(2e4),col="black",discrete=F,lwd=3)
#for (i in 1:1000) {
#  plot.cdf(rnorm(1e2),col=rgb(.5,.5,.5,.1),cex=.5,pch=20,add=T,print=F)
#  cat("\r Progress: ",i,"/",1000)
#}

ldir <- function(x,a,const=T) { # returns log dirichlet pdf
  a <- ifelse(a == 0, 1e-200, a)
  x <- ifelse(x == 0, 1e-200, x)
  num <- sum((a-1)*log(x))

  B <- 0
  if (const) {
    suma <- ifelse(sum(a) == 0, 1e-200, sum(a))
    B <- sum(lgamma(a)) - lgamma(suma)
  }
  out <- num - B

  out
}

gamma.mv2shsc <- function(m,v) {
  c("shape"=m^2/v, "scale"=v/m)
}

cdf2pdf <- function(x) {
  out <- c(x[1], x[-1] - x[-length(x)])
  ifelse(out == 0, 1e-200, out)
}
