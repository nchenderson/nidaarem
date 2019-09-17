UpdateStepLength <- function(par, X, y, lambda, stplngth) {
  beta.vec <- par

  stplngth <- stplngth/2
  X.beta <- X%*%beta.vec
  ## take step with "old" step-length
  beta.new <- SoftThresh(beta.vec + stplngth*crossprod(X, y - X.beta), lambda=lambda*stplngth)
  X.beta <- X%*%beta.new
  rr <- X.beta - 2*y
  f.old <- crossprod(X.beta, rr) + lambda*sum(abs(beta.new))
  beta.vec <- beta.new
  f.best <- f.old
  for(k in 1:3) {

    stp.new <- 2*stplngth
    beta.new <- SoftThresh(beta.vec + stp.new*crossprod(X, y - X.beta), lambda=lambda*stp.new)
    X.beta <- X%*%beta.new
    rr <- X.beta - 2*y
    f.new <- crossprod(X.beta, rr) + lambda*sum(abs(beta.new))
   # print(c(f.new,f.old))
    if(f.new > f.old) {
     # print("in break")
      ### monotonicity violation
      ### Use previous value of beta and step-length
      break
    }
    beta.vec <- beta.new

    f.old <- f.new
    stplngth <- stp.new
  }
  ans <- list(par=beta.vec, stplngth=stplngth)
  return(ans)
}


