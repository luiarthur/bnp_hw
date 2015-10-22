# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")

# Posterior Simulation ################
# y_i | G ~ G
# G ~ DP(a,G_0)

# True distributions for data

pD1 <- function(x) pnorm(x) # First True distribution (CDF) of y_i
pD2 <- function(x) { # 2. Second True distribution (CDF) of y_i
  .5*pnorm(x,-2.5,.5) + 0.3*pnorm(x,.5,.7) + 0.2*pnorm(x,1.5,2)
}

rD1 <- function(n) rnorm(n)
rD2 <- function (n) {
  r <- rmultinom(1:3,n,prob=c(.5,.3,.2))
  c(rnorm(r[1],-2.5,.5), rnorm(r[2],.5,.7), rnorm(r[3],1.5,2))
}
pData <- list(pD1,pD2)
rData <- list(rD1,rD2)

data.distribution <- list("cdf"=pData,"sampler"=rData)

# ams is alpha, m, and s. I'll have 8 rows. I want every 
# combination of big and small (or near and far).
sb <- rbind("a"= c(1,10), "m"= c(0,3), "s"= c(1,3))
ams <- matrix(c(sb[1,1],sb[2,1],sb[3,1],
                #sb[1,1],sb[2,1],sb[3,2],
                #sb[1,1],sb[2,2],sb[3,1],
                sb[1,1],sb[2,2],sb[3,2],
                sb[1,2],sb[2,1],sb[3,1],
                #sb[1,2],sb[2,1],sb[3,2],
                #sb[1,2],sb[2,2],sb[3,1],
                sb[1,2],sb[2,2],sb[3,2]), ncol=3,byrow=T)
colnames(ams) <- c("a","m","s")
ns <- c(20,200,2000)

# Fix:
# G0 = N(m,s)
# Pick one set: (m, s, a)


B <- 1000; J <- 1000
pG0.ms <- function(x,m,s) pnorm(x,m,s)
for (mod.num in 1:length(data.distribution[[1]])) { #2
  pD <- data.distribution$cdf[[mod.num]] 
  rD <- data.distribution$sampler[[mod.num]] 

  pdf(paste0("pdfs/dppost",mod.num,".pdf"),width=11,height=11)
  par(mfrow=c(length(ns),nrow(ams)))

  for (n.num in 1:length(ns)) { #3
    n <- ns[n.num]
    y <- rD(n)

    for (ams.num in 1:nrow(ams)){ #8
      a <- ams[ams.num,1] 
      m <- ams[ams.num,2] 
      s <- ams[ams.num,3] 
     
      pG0 <- function(x) pG0.ms(x,m,s)
      pG0_new <- function(x) {
        s <- sapply(x,function(w) sum(w >= y))
        (a * pG0(x) +  s) / (a + n)
      }

      # Draw from G | y
      xlim <- c(-9,9)
      G <- dp(N=B,pG0_new,a=a+n,xlim=xlim,J=J)
      EG <- apply(G$G,2,mean)

      dp.post.ci(G,ylab="Fn(y)",xlab="y",cex.main=1,lwd.EG=2,
                 main=bquote("("~alpha~"="~.(a)~", m ="~.(m)~", s="~.(s)~")"))
      plot.cdf(y,lwd=1,add=T,type="s")
      curve(pG0,add=T,col="red",lwd=1)  # prior G0
      curve(pD,add=T,col="green",lwd=1) # True Distribution of Data
      
      #legend("bottomright",lwd=2,bg=rgb(.9,.9,.9,.5),box.col=rgb(.9,.9,.9,.5),
      #       legend=c('Truth','Data',paste0('G0 = N(',m,',',s,')'),
      #                'Posterior:\n 95% C.I.','E[G|y]',''),cex=.5,
      #       col=c('green','black','red',rgb(.3,.3,.3),'blue','transparent'))
 

      cat("\r Model:", mod.num, "/2",
          " size: ",n.num, "/3",
          " ams:", ams.num, "/",nrow(ams))
    }
  }

  par(mfrow=c(1,1))
  dev.off()
}

