## This file contains basic functions that are useful elsewhere in the package
## Functions included are: the objective function for penalized logistic regression,
## objective functions for lasso, the lasso gradient descent updating function
## and the expit function

LogisticObjFn <- function(par, X, Xty, lambda) {
  X.beta <- as.vector(X%*%par)
  p1 <- sum(par*Xty)
  
  ## Do log-sum-exp trick for p2
  tmp <- X.beta > 0
  p2.pos <- sum(X.beta[tmp] + log1p(exp(-X.beta[tmp])))
  p2 <- sum(log1p(exp(X.beta[!tmp]))) + p2.pos
 
  p3 <- lambda*sum(abs(par))
  ans <- p1 - p2 - p3
  return(ans)
}

LassoObjFn_ngtp <- function(par, XtX, Xty, yty, lambda) {
  XX.beta <- as.vector(XtX%*%par)
  ans <- sum(par*Xty) - yty/2 - sum(par*XX.beta)/2 - lambda*sum(abs(par))
  return(ans)
}

LassoObjFn_pgtn <- function(par, X, Xty, yty, lambda) {
  X.beta <- as.vector(X%*%par)
  ans <- sum(par*Xty) - yty/2 - sum(X.beta*X.beta)/2 - lambda*sum(abs(par))
  return(ans)
}

GradDescFn <- function(par, X, y, lambda, stplngth) {
  X.beta <- X%*%par
  beta.vec <- SoftThresh(par + stplngth*crossprod(X, y - X.beta), lambda=lambda*stplngth)
  return(beta.vec)
}

expit <- function(x) {
  z <- x
  z[x > 0] <- 1/(1 + exp(-x[x>0]))
  z[x <= 0] <- exp(x[x<=0])/(1 + exp(x[x<=0]))
  return(z)
}

