system("mkdir -p ../../latex/img")
library(corrplot)
library(DPpackage) # For number 2
system("julia dpmm1.jl")
source("../../../R_Functions/plotPost.R",chdir=T)
source("../../../01_dp_priors/code/dp.R")
y <- read.table("../../dat/hw2.dat")[[1]]
alpha <- read.table("temp/out_alpha.dat")[[1]]
phi <- read.table("temp/out_phi.dat")[[1]]
mu <- read.table("temp/out_mu.dat")[[1]]
t2 <- read.table("temp/out_t2.dat")[[1]]
theta <- read.table("temp/out_theta.dat")

n <- length(y)
B <- nrow(theta)
burn <- round(B * .3)
ttheta <- tail(theta, B-burn)
talpha <- tail(alpha, B-burn)
tphi <- tail(phi, B-burn)
tmu <- tail(mu, B-burn)
tt2 <- tail(t2, B-burn)

#mtheta <- apply(tail(theta,B-burn),2,mean)
#plot(y,pch=20,col="grey",cex=2)
#points(rnorm(n,mtheta,sqrt(mphi)),pch=20,col='blue')


# Posterior for individual Parameters:
#pdf("../../latex/img/posts1alphaorig.pdf",width=10,height=10)
par(mfrow=c(2,2))
  plot.post(talpha,main=bquote(alpha),cex.main=2)
  plot.post(tphi,main=bquote(phi),cex.main=2)
  plot.post(tmu,main=bquote(mu),cex.main=2)
  plot.post(tt2,main=bquote(tau^2),cex.main=2)
par(mfrow=c(1,1))
#dev.off()

# Number of Clusters:
clus <- apply(ttheta,1, function(x) length( unique (x) ) )
tabclus <- table(clus)

#pdf("../../latex/img/nstar1big.pdf")
plot(tabclus / sum(tabclus),type='h',bty='n',fg='grey',xlab='Clusters',
     ylab = 'Proportions',lwd=3,col="blue")
#dev.off()
mean(clus) # E[n*|...]
var(clus)

# Trace Plot for Thetas ####################################

plot(ttheta[,1],pch=20,fg='grey',bty="n",cex=.1,col=rgb(.4,.4,.4,.1),
     ylim=range(ttheta), main=bquote("Trace Plot of "~theta))
for ( k in 2:n ) {
  points(ttheta[,k],col=rgb(.4,.4,.4,.1),pch=20,fg='grey',bty="n",cex=.1)
}

# Cluster Behavior: ###################################
mtheta <- apply(ttheta,2,mean)
qtheta <- apply(theta,2,function(x)quantile(x,c(.025,.5,.957)))
order_y <- order(mtheta)
Clus <- matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n) {
    Clus[i,j] <- sum(ttheta[,i] == ttheta[,j]) / (B-burn)
  }
}
w2b <- rgb(seq(1,0,len=100),seq(1,0,len=100),1)
image(Clus[order_y,order_y],xaxt="n",yaxt="n",
      col=w2b)
axis(1,tck=0,at=1:n/n,lab=order_y,cex.axis=.4,las=2)
axis(2,tck=0,at=1:n/n,lab=order_y,cex.axis=.4,las=2)
corrplot(Clus,order="hclust",method='square',insig='blank')

clust <- hclust(dist(mtheta))
image(Clus[clust$order,clust$order],col=w2b)
dev.off()

ord <- order(y)

#pdf("../../latex/img/clus1.pdf")
plot(y[ord],type='p',col='grey30',cex=.5,fg='grey',bty='l',ylab="y",xlab="Ordered Observations",pch=20)
for (i in 1:n)
  segments(i,qtheta[1,ord[i]],i,qtheta[3,ord[i]],col='pink')
points(qtheta[2,ord],cex=.3,col='purple')
legend("topleft",legend=c("Data","Posterior Median","2.5 - 97.5% Quantiles"),
       col=c("grey30","purple","pink"),lwd=3,bty='n')
#dev.off()

# Posterior Predictive: ###############################
prior.a <- rgamma(10000,2,4) #2,4
prior.phi <- 1 / rgamma(10000,3,3) # 3,15
prior.mu <- rnorm(10000,0,3) #0,3
prior.t2 <- 1 / rgamma(10000,3,15) #3,6
y.prior <- apply(matrix(1:10000),1,function(b) {
  cat("\r Progress: ",b)
  th <- dp(1,prior.a[b],function(x) pnorm(x,prior.mu[b],sqrt(prior.t2[b])),
           xlim=c(-10,10))
  ind <- which(th$G > runif(1))[1]
  th <- th$x[ind]
  rnorm(1,th,sqrt(prior.phi[b]))
})

y.pred <- apply(matrix(1:(B-burn)),1,function(b) {
            th <- 0

            if (talpha[b] / (talpha[b]+n) > runif(1)) {
              th <- rnorm(1,tmu[b],sqrt(tt2[b]))
            } else {
              th <- sample(ttheta[b,],1)
            }

            cat("\r Progress: ",b)
            rnorm(1, as.numeric(th), sqrt(tphi[b]))
      })
f <- function(x) {.2*dnorm(x,-5,1)+.5*dnorm(x) + .3*dnorm(x,3.5,1)}
# Check soultion with DPpackage: ######################################
#prior1 <- list(a0=1,b0=1,m2=0,s2=1,psiinv1=diag(.5,1),nu1=1,tau1=1,tau2=1)
#dp.den <- DPdensity(y,prior=prior1,
#                    mcmc=list(nburn=burn,nsave=B-burn,ndisplay=B*.01,nskip=0),
#                    state=NULL,status=T)
#dims <- dim(dp.den$save.state$randsave)
##Posterior Predictive:
#lines(density(dp.den$save.state$randsave[,dims[2]]),col='red',lwd=3)

# 2: Model with DPpackage Problem ####################
prior2 <- list(a0=3,b0=3,m2=0,s2=3,psiinv2=solve(diag(.5,1)),nu1=1,nu2=2,tau1=6,tau2=6)
dp.den <- DPdensity(y,prior=prior2,
                    mcmc=list(nburn=burn,nsave=B-burn,ndisplay=B*.01,nskip=0),
                    state=NULL,status=T)
ind <- which(colnames((dp.den$save.state$randsave)) == "y (Prediction)")
# Posterior Predictive:
#pdf("../../latex/img/postpred1HighT2.pdf")
plot(density(y.pred),col='blue',ylim=c(0,.2),lwd=2,bty='n',fg='grey',
     main="Posterior Predictive")
lines(density(y),col='grey',lwd=3)
curve(f, add=T,col='green',lwd=2)
lines(density(y.prior),lwd=2,col="gold")
#lines(density(dp.den$save.state$randsave[,ind]),col='red',lwd=3)
legend("topleft",legend=c("Posterior Predictive","Data","Truth","Prior Predictive"),
       col=c("blue","grey","green","gold"),bty='n',lwd=2)
#dev.off()



