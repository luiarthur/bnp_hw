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

# MDP (Ferguson)
#N = 2
#rA <- r.a.prior
#pG <- function(x) pnorm(x)
#xlim = c(-3,3)
#n <- 100
mdp <- function(N=1,rA,pG,xlim=c(0,1),n=100) {
  x <- seq(xlim[1],xlim[2],length=n)
  x <- c(-1e1000,x)
  a <- rA(N) # prior for alpha
  adG0 <- apply(matrix(a),1, function(aa) aa*pG(x[-1]) - aa*pG(x[-length(x)]))
  #out <- rdir(N,c(adG0))

  y <- matrix(rgamma(n*N, adG0, 1), N, n, byrow=T)
  rowsums <- y %*% rep(1,n)
  out <- y / as.vector(rowsums)

  G <- t(apply(out,1,cumsum))
  list("G"=G, "x"= x[-1],"a"=a)
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
}

dp.post.extra <- function(X,col.lines=rgb(.4,.4,.4,.1),xlim.def=range(X$x),...) {
  dp.post(X,...)
  xx <- X$x

  lines(xx,apply(X$G,2,var),col="red",lwd=2)

  legend("right",legend=c("Random draws from the DP",
                          "A particular draw from the DP",
                          "E[G]","Var[G]"),
         col=c("grey","black","blue","red"),lwd=3,
         bg=rgb(.9,.9,.9,.5),box.col=rgb(.9,.9,.9,.5) )

  minor <- function() {
    plot(xx,apply(X$G,2,var),col="red",lwd=1,bty="n",type="l",
         col.axis=rgb(.3,.3,.3),fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),
         col.main=rgb(.3,.3,.5),cex.axis=.5)
  }
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

# MAIN:


# a) DP Sethuraman & Ferguson

# Setharuman
pdf("pdfs/sethDP.pdf",width=15,height=9)
par(mfrow=c(1,3))
avec <- c(1,10,100)
for (a in avec) {
  gs <- dp_stickbreak(N=1000, a=a,rG=function(n) rnorm(n), xlim=c(-3,3), eps=1e-4)
  #main <- bquote(paste("DP(",.(al),"G"[0],")"))
  dp.post.extra(gs,col.lines=rgb(.4,.4,.4,.05),ylab="F(x)",xlab="x",
          main=bquote("G ~ DP("~.(a)~","~G[0]~")"~" - Sethuraman's Construction"))
}
par(mfrow=c(1,1))
dev.off()

# Ferguson
pdf("pdfs/fergusonDP.pdf",width=15,height=9)
par(mfrow=c(1,3))
for (a in avec) {
  gf <- dp(N=1000, a=a,pG=function(n) pnorm(n), xlim=c(-3,3))
  dp.post.extra(gf,col.lines=rgb(.4,.4,.4,.05),ylab="F(x)",xlab="x",
          main=bquote("G ~ DP("~.(a)~","~G[0]~")"~" - Ferguson's Construction"))
}
par(mfrow=c(1,1))
dev.off()

# b) Mixture of DPs prior (MDP):  
# G|a ∼ DP[a, N(0,1)]
#   a ~ .3G(3,2) + .6G(9,6) + .1G(1,9)

mdp.post <- function(X,col.lines=rgb(.4,.4,.4,.1),xlim.def=range(X$x),...) {
  dp.post(X,...)
  xx <- X$x
  lines(xx,apply(X$G,2,var),col="red",lwd=2)

  legend(0,.2,legend=c("Random draws from the DP",
                          "A particular draw from the DP",
                          "E[G]","Var[G]",expression(paste("Density of ",alpha))),
         col=c("grey","black","blue","red","orange"),lwd=3,
         bg=rgb(.9,.9,.9,.5),box.col=rgb(.9,.9,.9,.5) )

  if (!is.null(gmdp$a)) {
    minor <- function() {
      plot(density(gmdp$a),col="orange",lwd=1,bty="n",type="l",main="",
           col.axis=rgb(.3,.3,.3),fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),
           col.main=rgb(.3,.3,.5),cex.axis=.5)
    }
    plot.in.plot(minor,"topleft")
  }
}

r.a.prior <- function(n) {
  count.of.each.dist <- c(rmultinom(1,n,c(.3,.6,.1)))
  r1 <- rgamma(count.of.each.dist[1],3,scale=2)
  r2 <- rgamma(count.of.each.dist[2],9,scale=6)
  r3 <- rgamma(count.of.each.dist[3],1,scale=9)
  out <- c(r1,r2,r3)
  out
}

rA <- list(function(n) rgamma(n,3,sc=2),
           function(n) rgamma(n,6,sc=4),
           function(n) rgamma(n,10,sc=1))
pdf("pdfs/priorMDP.pdf",width=15,height=9)
par(mfrow=c(1,3))
for (ra in rA) {
  gmdp <- mdp(N=1000, rA=ra, pG=function(n) pnorm(n), xlim=c(-3,3))
  mdp.post(gmdp,main=bquote("G|a ~ DP(a,"~G[0]~")"))
}
par(mfrow=c(1,1))
dev.off()
