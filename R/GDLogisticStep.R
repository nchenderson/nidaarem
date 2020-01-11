GDLogisticStep <- function(par, X, Xty, lambda, stplngth) {
  phat <- expit(X%*%par)
  xnew <- SoftThresh(par + stplngth*(Xty - crossprod(X, phat)), lambda=lambda*stplngth)
  return(xnew)
}