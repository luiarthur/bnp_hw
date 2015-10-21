# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/plotPost.R",chdir=T)
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")
set.seed(1)

pD1 <- function(x) ppois(x,5) # First True distribution (CDF) of y_i
pD2 <- function(x) .7*ppois(x,3) + 0.3*ppois(x,11)

rD1 <- function(n) rpois(n,5)
rD2 <- function (n) {
  r <- rmultinom(1:2,n,prob=c(.7,.3))
  c(rpois(r[1],3),rpois(r[2],11))
}

csa.12 <- c(6,6)
csl.12 <- c(2,2)
pData <- list(pD1,pD2)
rData <- list(rD1,rD2)

data.distribution <- list("cdf"=pData,"sampler"=rData)

n <- 300
B <- 5000
burn <- round(B * .2)
xlim <- 0:30

# Start here
for (mod.num in 1:length(data.distribution[[1]])){
  y <- rData[[mod.num]](n)
  n.uniq <- length(unique(y))

  table.y <- table(y)
  y.star <- as.numeric(names(table.y))
  nj <- as.numeric(table.y)

  csa <- csa.12[mod.num] # cand sig for alpha
  csl <- csl.12[mod.num] # cand sig for lambda
  
  ppa <- gamma.mv2shsc(5,5) # Prior for alpha is Gamma(5,1) because no idea
  ppl <- gamma.mv2shsc(5,5) # Prior for lambda is Gamma(5,1) because center is 5
  lpa <- function(x) dgamma(x,ppa[1],sc=ppa[2],log=T) #log prior for alpha
  lpl <- function(x) dgamma(x,ppa[1],sc=ppa[2],log=T) #log prior for lambda
  
  K <- length(xlim)
  G <- matrix(0,B,K)
  lg0 <- function(x,l) dpois(x,l,log=T)
  pG0 <- function(x,l) ppois(x,l)
  dG0 <- function(x,l) dpois(x,l)
  pGn <- pG0 # this is an initialization
  G[1,] <- dp(N=1,a=1, pG=function(x) pG0(x,1), xlim=xlim)$G
  lam <- rep(1,B)
  a <- rep(1,B)
  acc.l <- acc.a <- 0
  
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
    canda <- rnorm(1,a[b],csa) # Normal
    if (canda > 0) {
      lla <- function(aa) {
        n.uniq * log(aa) - (lgamma(aa+n) - lgamma(aa)) +
        sum( lgamma(aa*exp(lg0(y.star,lam[b])) + nj) - 
             lgamma(aa*exp(lg0(y.star,lam[b])) + 1))
      }

      lga1 <- lla(canda) + lpa(canda)
      lga2 <- lla(a[b]) + lpa(a[b])
      lr <- lga1 - lga2
  
      if (lr > log(runif(1))) {
        a[b] <- canda
        acc.a <- acc.a + 1
      }
    }
  
    # Update lambda
    candl <- rnorm(1,lam[b],csl) # normal proposal
    if (candl > 0) {
      lll <- function(l) {
        sum( lg0(y.star,l) + 
             lgamma(a[b]*exp(lg0(y.star,l)) + nj) - 
             lgamma(a[b]*exp(lg0(y.star,l)) + 1 ))
      }

      lgl1 <- lll(candl)   + lpl(candl)
      lgl2 <- lll(lam[b]) + lpl(lam[b])
      lr <- lgl1 - lgl2
  
      if (lr > log(runif(1))) {
        lam[b] <- candl
        acc.l <- acc.l + 1
      }
    }
  
    cat("\r ", b / B)
  }
  #out <- list("G"=G, "a"=a, "lamda"=lam, "acc.a"=acc.a, "acc.l"=acc.l)
  
  #G
  pdf(paste0("pdfs/postMDP",mod.num,".pdf"),w=10,h=10)
  #layout(matrix(c(1,1,2,3,4,5),2,byrow=T))
  #layout(matrix(c(1,3,4,2),2,byrow=T))
  par(mfrow=c(2,2))
  dp.post.ci(list("G"=G[-c(1:burn),],"x"=xlim),ylab="Fn(y)",
             xlab="y",type.EG="p",pch=20,cex.EG=3,
             main=bquote(" G | y,"~alpha~","~lambda ),
             cex.main=2,EG.col=rgb(.3,.3,.7),dens=T,
             ylim.def=c(0,.3),lwd.ci=3)
  pmf <- plot.pmf(y,type="p",add=T,pch=20,col="green",cex=1.5)
  legend("topright",
         col=c("green",rgb(.3,.3,.7),"red","grey"),text.col=rgb(.3,.3,.4),
         legend=c("Data","E[G|y]",expression(paste("[","y"^"new","|y]")),"95% C.I."),lwd=3,
         bg=rgb(.9,.9,.9,.5),box.col=rgb(.9,.9,.9,.5),
         cex=1.5)

  # Posterior Predictive
  postpred <- apply(matrix(tail(1:B,B-burn)),1,function(i) {
                    #prob <- rdir( 1, a[i] * dG0(xlim,lam[i]) )
                    dG <- c(G[i,1],G[i,-1]-G[i,-length(xlim)]) 
                    prob <- rdir( 1, a[i] * dG )
                    sample(xlim,1,prob=prob)
         })
  pmf.postpred <- plot.pmf(postpred,type="p",add=T,pch=20,col=rgb(1,0,0,.5),cex=2)

  # alpha posterior
  plot.post(tail(a,B-burn),
            main=bquote("Posterior & Trace for"~alpha~" ("~.(round(100*acc.a/(B)))~"% acceptance)"))

  # lambda posterior
  plot.post(tail(lam,B-burn),
            main=bquote("Posterior & Trace for"~lambda~" ("~.(round(100*acc.l/(B)))~"% acceptance)"))

  # alpha lambda joint
  plot(tail(lam,B-burn),tail(a,B-burn),type="p",main=bquote(alpha~" vs "~lambda), col=rgb(.2,.2,((burn+1):B)/B,.2),cex=.4,bty="n",
       ylab=bquote(alpha),xlab=bquote(lambda),fg='grey')

  dev.off()

}
#system("cd ../latex; comptex")
#source("posteriorMDP.R")
