y <- read.table("../../dat/hw2.dat")[[1]]
alpha <- read.table("temp/out_alpha.dat")
phi <- read.table("temp/out_phi.dat")
mu <- read.table("temp/out_mu.dat")
t2 <- read.table("temp/out_t2.dat")
eta <- read.table("temp/out_eta.dat")
theta <- read.table("temp/out_theta.dat")

mtheta <- apply(theta,2,mean)
plot(y,pch=20,col="grey",cex=2)
points(mtheta,pch=20,col='blue')
