auxGibbsR <- function(y, a=1, s=1, cs =3, B=10000) {
  n <- length(y)
  theta <- matrix(0,B,n)
  rG0 <- rnorm(B,0,3)
  p <- 0

  for (b in 2:B) {
    theta[b,] <- theta[b-1,]
    for (i in 1:n) {
      # Update theta_i | theta_{-1}
      tb <- theta[b,]
      ti <- tb[i]
      ut_xi <- unique( tb[-i] )
      k <- length( ut_xi )

      if ( ti %in% ut_xi ) p <- rG0[b]
      probs <- c( sum(ut_xi == tb) * dnorm(y[i],ut_xi) , a * dnorm(y[i],p) )

      theta[b,i] <- sample( c(ut_xi,p), 1, prob=probs)
    }

    # Update theta | y
    tb <- theta[b,]
    ut <- unique(tb)
    for (u in ut) {
      ind <- which(tb == u)
      cand <- rnorm(1,u,cs)
      lg1 <- sum(dnorm(y[ind],cand,log=T)) + lg0(cand)
      lg2 <- sum(dnorm(y[ind],u,log=T)) + lg0(u)
      lg <- lg1 - lg2

      if (lg > log(runif(1)) ) {
        #acc_t <- acc_t + 1
        theta[b,i] <- cand
      }
    }
    # Print Progress:
    if (b %% (B/20)==0) cat("\rProgress: ", round(100*b/B), "%")
  } # end of R_gibbs

  theta
}
