# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")

# Posterior Simulation ################
# y_i | G ~ G
# G ~ DP(a,G_0)

pT1 <- function(x) pnorm(x) # First True distribution (CDF) of y_i
pT2 <- function(x) { # 2. Second True distribution (CDF) of y_i
  .5*pnorm(x,-2.5,.5) + 0.3*pnorm(x,.5,.7) + 0.2*pnorm(x,1.5,2)
}

# Fix:
# G0 = N(m,s)
# Pick one set: (m, s, a)
a <- 5 # Later, put a prior. Large => More faith in prior.
m <- 7
s <- 3
pG0 <- function(x) pnorm(x,m,s)

# Select sample size:
i <- 2 # 1,2,3
n <- c(20,200,2000)
y1 <- rnorm(n[i])
r <- rmultinom(1:3,n[i],prob=c(.5,.3,.2))
y2 <- c( rnorm(r[1],-2.5,.5), rnorm(r[2],.5,.7), rnorm(r[3],1.5,2))

# Gibbs
B <- 1000 # Big number for MCMC
J <- 1000 # Number of obs for each DP draw, should be LARGE

y <- y2 # change this y1 or y2
pG_new <- function(x) {
  s <- sapply(x,function(w) sum(w >= y))
  (a * pG0(x) +  s) / (a + n[i])
}

# Draw from G | y
xlim <- c(-5,10)
G <- dp(N=B,pG_new,a=a+n[i],xlim=xlim,J=J)
EG <- apply(G$G,2,mean)

  plot(0,cex=0,ylim=c(0,1),xlim=xlim,
       #main=bquote(alpha~"="~.(a)~","~G[0]~"="~"N("~.(m)~","~.(s)~");"~" n ="~.(n[i])),
       main=bquote("G ~ DP("~.(a)~","~"N("~.(m)~','~.(s)~"));"~" n ="~.(n[i])),
       ylab="Fn(y)",xlab="y", bty="n",las=1, col.axis=rgb(.3,.3,.3),
       fg=rgb(.8,.8,.8),col.lab=rgb(.3,.3,.5),col.main=rgb(.3,.3,.4))
  lines(ecdf(y),cex=.5)
  curve(pG0,add=T,col="red",lwd=3)
  lines(G$x,EG,col=rgb(0,0,1,.3),lwd=10)
  curve(pT2,add=T,col="green",lwd=3) # change pT1, pT2
  qG <- apply(G$G,2,function(x) quantile(x,c(.025,.5,.975)))
  glo <- qG[1,]
  ghi <- qG[3,]
  color.btwn(G$x,glo,ghi,-100,100,col.area=rgb(.2,.2,.2,.5))

  legend("bottomright",lwd=3,bg=rgb(.9,.9,.9,.5),box.col=rgb(.9,.9,.9,.5),
         legend=c('Truth','Data',paste0('G0 = N(',m,',',s,')'),
                  'Posterior:\n 95% C.I.','E[G|y]',''),
         col=c('green','black','red',rgb(.3,.3,.3),'blue','transparent'))

