# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")

pD1 <- function(x) ppois(x,5) # First True distribution (CDF) of y_i
pD2 <- function(x) .7*ppois(x,3) + 0.3*ppois(x,11)

rD1 <- function(n) rpois(n,5)
rD2 <- function (n) {
  r <- rmultinom(1:2,n,prob=c(.7,.3))
  c(rpois(r[1],3),rpois(r[2],11))
}
pData <- list(pD1,pD2)
rData <- list(rD1,rD2)

data.distribution <- list("cdf"=pData,"sampler"=rData)

n <- 300
B <- 10000
xlim <- 0:30

# Start here
y <- rD1(n) 
csa <- 1 # cand sig for alpha
csl <- 1 # cand sig for lambda

lpa <- function(x) dgamma(x,1,1,log=T) #log prior for alpha
lpl <- function(x) dgamma(x,1,1,log=T) #log prior for lambda
lq <- function(x,m,v) {
  ss <- gamma.mv2shsc(m,v)
  dgamma(x,ss[1],sc=ss[2],log=T)
}

K <- length(xlim)
G <- matrix(0,B,K)
pG0 <- function(x,lam) ppois(x,lam)
pGn <- pG0 # this is an initialization
G[1,] <- dp(N=1,a=1, pG=function(x) pG0(x,1), xlim=xlim)$G
lam <- rep(1,B)
a <- rep(1,B)
acc.l <- acc.a <- 0

# TEST
#xlim <- 0:20
#xlim <- c(0,20)
#M <- dp(N=1000, a=2, pG=function(x) ppois(x,3), xlim=xlim)
#M <- dp(N=1000, a=3+n, pG=function(x) pGn(x,3), xlim=xlim)
#dp.post(M)
#lines(ecdf(y),lwd=3)
#lines(M$x,apply(M$G,2,function(x) mean(x,na.rm=T)),col="blue",type="s",lwd=3)
#lines(M$x, ppois(M$x,3), col="green",lwd=3,type="s")

for (b in 2:B) {
  lam[b] <- lam[b-1]
  a[b] <- a[b-1]
  G[b,] <- G[b-1,]

  # Update G
  pGn <- function(x,lamb) {
    s <- sapply(x,function(w) sum(w >= y))
    (a[b] * pG0(x,lamb) +  s) / (a[b] + n)
  }
  G[b,] <- dp(N=1, a=a[b]+n, pG=function(x) pGn(x,lam[b]), xlim=xlim)$G
  
  # Update alpha
  #ss <- gamma.mv2shsc(a[b],csa) # gamma
  #cand <- rgamma(1,ss[1],sc=ss[2]) # gamma proposal
  cand <- rnorm(1,a[b],csa) # Normal
  if (cand > 0) {
    pg <- cdf2pdf(G[b,])
    pg01 <- c(pG0(xlim[1],lam[b]),pG0(xlim[-1],lam[b])-pG0(xlim[-length(xlim)],lam[b]))
    pg02 <- c(pG0(xlim[1],lam[b]),pG0(xlim[-1],lam[b])-pG0(xlim[-length(xlim)],lam[b]))

    lga1 <- ldir(pg,cand*pg01) + lpa(cand)
    lga2 <- ldir(pg,a[b]*pg02) + lpa(a[b])
    lqa1 <- 0 #lq(cand,a[b],csa) #0 Normal
    lqa2 <- 0 #lq(a[b],cand,csa) #0 Normal
    lr <- lga1 - lqa1 - lga2 + lqa2

    if (lr > log(runif(1))) {
      a[b] <- cand
      acc.a <- acc.a + 1
    }
  }

  # Update lambda
  cand <- rnorm(1,lam[b],csl) # normal proposal
  #ss <- gamma.mv2shsc(lam[b],csl) # gamma proposal
  #cand <- rgamma(1,ss[1],sc=ss[2]) # gamma proposal
  if (cand > 0) {
    pg01 <- c(pG0(xlim[1],cand),   pG0(xlim[-1],cand)   - pG0(xlim[-length(xlim)],cand))
    pg02 <- c(pG0(xlim[1],lam[b]), pG0(xlim[-1],lam[b]) - pG0(xlim[-length(xlim)],lam[b]))

    lgl1 <- ldir(pg,a[b]*pg01) + lpl(cand)
    lgl2 <- ldir(pg,a[b]*pg02) + lpl(lam[b])
    lql1 <- 0 #lq(cand,lam[b],csl) #0
    lql2 <- 0 #lq(lam[b],cand,csl) #0
    lr <- lgl1 - lql1 - lgl2 + lql2

    if (lr > log(runif(1))) {
      lam[b] <- cand
      acc.l <- acc.l + 1
    }
  }

  cat("\r ", b / B)
}
#list("G"=G, "a"=a, "lamda"=lam, "acc.a"=acc.a, "acc.l"=acc.l)

#G
layout(matrix(c(1,1,1,2,3,4),2,byrow=T))
dp.post.ci(list("G"=G[-c(1:8500),],"x"=xlim),ylab="Fn(y)",xlab="y",cex.main=.9,type.EG="p",pch=20,cex.EG=3)
plot.cdf(y,type="s",add=T)
plot(a[-c(1:200)],type="l",main=paste("alpha", round(100*acc.a/(B),5),"%"))
plot(a[-c(1:200)],lam[-c(1:200)],type="p",main="alpha vs. lambda", col=rgb(.2,.2,(201:B)/B,.2),cex=.4)
plot(lam[-c(1:200)],type="l",main=paste("lambda", round(100*acc.l/(B),5),"%"))

#source("posteriorMDP.R")
