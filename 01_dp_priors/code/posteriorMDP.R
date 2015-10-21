# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/plotPost.R",chdir=T)
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")
#set.seed(1)

pD1 <- function(x) ppois(x,5) # First True distribution (CDF) of y_i
pD2 <- function(x) .7*ppois(x,3) + 0.3*ppois(x,11)

rD1 <- function(n) rpois(n,5)
rD2 <- function (n) {
  r <- rmultinom(1:2,n,prob=c(.7,.3))
  c(rpois(r[1],3),rpois(r[2],11))
}

csa.12 <- c(5,5)
csl.12 <- c(.5,1)
pData <- list(pD1,pD2)
rData <- list(rD1,rD2)

data.distribution <- list("cdf"=pData,"sampler"=rData)

n <- 300
B <- 15000
burn <- B * .2
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
  
  ppa <- gamma.mv2shsc(5,5)
  ppl <- gamma.mv2shsc(5,5)
  lpa <- function(x) dgamma(x,ppa[1],sc=ppa[2],log=T) #log prior for alpha
  lpl <- function(x) dgamma(x,ppa[1],sc=ppa[2],log=T) #log prior for lambda
  lq <- function(x,m,v) {
    ss <- gamma.mv2shsc(m,v)
    dgamma(x,ss[1],sc=ss[2],log=T)
  }
  
  K <- length(xlim)
  G <- matrix(0,B,K)
  lg0 <- function(x,lam) dpois(x,lam,log=T)
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
      #pg <- cdf2pdf(G[b,])
      #pg01 <- c(pG0(xlim[1],lam[b]),pG0(xlim[-1],lam[b])-pG0(xlim[-length(xlim)],lam[b]))
      #pg02 <- c(pG0(xlim[1],lam[b]),pG0(xlim[-1],lam[b])-pG0(xlim[-length(xlim)],lam[b]))
      #lga1 <- ldir(pg,cand*pg01) + lpa(cand)
      #lga2 <- ldir(pg,a[b]*pg02) + lpa(a[b])
      
      lla <- function(aa) {
        n.uniq * log(aa) - lgamma(aa+n) -lgamma(aa) +
        sum( (nj-1) * log(aa*exp(lg0(y.star,lam[b])) + 1) )
      }

      lga1 <- lla(cand) + lpa(cand)
      lga2 <- lla(a[b]) + lpa(a[b])

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
      #pg01 <- c(pG0(xlim[1],cand),   pG0(xlim[-1],cand)   - pG0(xlim[-length(xlim)],cand))
      #pg02 <- c(pG0(xlim[1],lam[b]), pG0(xlim[-1],lam[b]) - pG0(xlim[-length(xlim)],lam[b]))
  

      lll <- function(l) {
        sum( lg0(y.star,l) + 
            (nj-1) * log(a[b]*exp(lg0(y.star,l)) + 1) )
      }

      lgl1 <- lll(cand)   + lpl(cand)
      lgl2 <- lll(lam[b]) + lpl(lam[b])
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
  postpred <- apply(matrix(tail(1:n,B-burn)),1,function(i) {
                    dG0 <- c(pG0(xlim[1],lam[i]), pG0(xlim[-1],lam[i]) - pG0(xlim[-length(xlim)],lam[i]))
                    prob <- rdir(1,a[i]*dG0)
                    sample(xlim,1,prob=prob)
         })
  pmf.postpred <- plot.pmf(postpred,type="p",add=T,pch=20,col=rgb(1,0,0,.5),cex=2)

  # alpha posterior
  plot.post(tail(a,B-burn),
            main=bquote("Posterior & Trace for"~alpha~" ("~.(round(100*acc.a/(B),5))~"% acceptance)"))

  # lambda posterior
  plot.post(tail(lam,B-burn),
            main=bquote("Posterior & Trace for"~lambda~" ("~.(round(100*acc.l/(B),5))~"% acceptance)"))

  # alpha lambda joint
  plot(a[-c(1:burn)],lam[-c(1:burn)],type="p",main=bquote(alpha~" vs "~lambda), col=rgb(.2,.2,((burn+1):B)/B,.2),cex=.4,bty="n",
       ylab=bquote(alpha),xlab=bquote(lambda),fg='grey')

  dev.off()

}
#system("cd ../latex; comptex")
#source("posteriorMDP.R")
