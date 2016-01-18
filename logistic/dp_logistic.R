# Method 1:
#  mu = (mu_z, mu_x)
#   S = S_{zz}, S_{zx}
#       S_{xz}, S_{xx}

# [y_i , x_i1, x_ik] ~ Bern(y_i, p_i) x N_k(mu_x, S_{xx}), 
#   for S_{zx} = 0, 
#   p_i = mu_z / sqrt(S_zz), where S_zz is a constant 
# mu_x 

# Method 2: (Easier)
# yi | pi ~ Bern(pi)
#      pi  = (1+exp(-b0_i - xib)

# Generate Data:
set.seed(1)
n <- 1000
f <- function(x) dnorm(x,0,1) + dnorm(x,3,1) + .3
x <- runif(n,-5,10)
fx <- f(x)
y <- rbinom(n,1,fx)

mod.glm <- glm(y~x, family="binomial")
pred.glm <- predict(mod.glm,type="response")

plot(x,pred.glm,pch=20,cex=.1,ylim=c(0,1),bty='n')
points(x,fx,pch=20,col='green',cex=.1)
legend("topright",legend=c("Data","GLM","BNP"),col=c('green','black',"blue"),lwd=3,bty='n')


