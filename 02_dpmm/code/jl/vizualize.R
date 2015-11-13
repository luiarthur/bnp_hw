library(corrplot)
system("julia dpmm1.jl")
source("../../../R_Functions/plotPost.R",chdir=T)
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
par(mfrow=c(2,2))
  plot.post(talpha,main="alpha")
  plot.post(tphi,main="phi")
  plot.post(tmu,main="mu")
  plot.post(tt2,main="t2")
par(mfrow=c(1,1))
dev.off()

# Number of Clusters:
clus <- apply(ttheta,1, function(x) length( unique (x) ) )
tabclus <- table(clus)
plot(tabclus / sum(tabclus),type='h',bty='n',fg='grey',xlab='Clusters',
     ylab = 'proportions',
     main=paste("mean =",round(mean(clus),1)),lwd=3,col="blue")

# Trace Plot for Thetas ####################################
plot(ttheta[,1],pch=20,fg='grey',bty="n",cex=.1,col=rgb(.1,.1,.1,.1),
     ylim=range(ttheta), main=bquote("Trace Plot of "~theta))
for ( k in 2:n ) {
  points(ttheta[,k],col=rgb(.1,.1,.1,.1),pch=20,fg='grey',bty="n",cex=.1)
}

# Cluster Behavior: ###################################
order_y <- order(y)
Clus <- matrix(0,n,n)
for (i in 1:n) {
  for (j in 1:n) {
    Clus[i,j] <- sum(ttheta[,i] == ttheta[,j])
  }
}
image(Clus[order_y,order_y],xaxt="n",yaxt="n",
      col=rgb(seq(1,0,len=100),seq(1,0,len=100),1))
axis(1,tck=0,at=1:n/n,lab=order_y,cex.axis=.4,las=2)
axis(2,tck=0,at=1:n/n,lab=order_y,cex.axis=.4,las=2)

# Posterior Predictive: ###############################

