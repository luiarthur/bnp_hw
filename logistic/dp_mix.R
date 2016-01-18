n <- 100

beta <- c(0,1,6,-1,-12,2,20,-2)

x1 <- runif(n,0,3)
x2 <- runif(n,3,6)
x3 <- runif(n,6,8)
x4 <- runif(n,8,10)

x.list <- list(x1,x2,x3,x4)
x <- unlist(x.list)

y <- rep(0,n*4)
for (i in 1:4) {
  for (j in 1:n) {
  y[(i-1)*n+j] <- x.list[[i]][j] * beta[2*i] + beta[2*i-1]
  }
}

y <- rnorm(4*n,y,1)

plot(x,y)



