# Posterior Inference for one-sample problems using DP Priors
system("mkdir -p pdfs") # mkdir pdfs if it doesn't exist.
source("../../R_Functions/plotinplot.R")
source("../../R_Functions/colorUnderCurve.R")
source("dp.R")

pD1 <- function(x) ppois(x,5) # First True distribution (CDF) of y_i
pD2 <- function(x) .7*ppois(x,3) + 0.3*ppois(x,11)

rD1 <- function(n) rpois(n,5)
rD2 <- function (n) {
  r <- rmultinom(1:2,n,prob=c(.7,.3))
  c(rpois(r[1],3),rpois(r[2],11))
}
pData <- list(pD1,pD2)
rData <- list(rD1,rD2)

data.distribution <- list("cdf"=pData,"sampler"=rData)

n <- 300

B <- 1000; J <- 1000

