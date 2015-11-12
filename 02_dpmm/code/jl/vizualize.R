y <- read.table("../../dat/hw2.dat")[[1]]
alpha <- read.table("temp/out_alpha.dat")[[1]]
phi <- read.table("temp/out_phi.dat")[[1]]
mu <- read.table("temp/out_mu.dat")[[1]]
t2 <- read.table("temp/out_t2.dat")[[1]]
eta <- read.table("temp/out_eta.dat")[[1]]
theta <- read.table("temp/out_theta.dat")

B <- nrow(theta)
burn <- B * .3
mtheta <- apply(tail(theta,B-burn),2,mean)
plot(y,pch=20,col="grey",cex=2)
points(mtheta,pch=20,col='blue')


plot(density(tail(alpha,B-burn)))
plot(density(tail(phi,B-burn)))
plot(density(tail(mu,B-burn)))
plot(density(tail(t2,B-burn)))
plot(alpha,type="l")
plot(phi,type="l")
plot(mu,type="l")
plot(t2,type="l")

clus <- apply(theta,1, function(x) length( unique (x) ) )
plot(table(clus),type='h')
mean(clus)
