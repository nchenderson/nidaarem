LassoObjFn <- function(bbeta, X, y, lambda, stplngth) {
  X.beta <- as.vector(X%*%bbeta)
  ans <- sum((y - X.beta)*(y - X.beta))/2 + lambda*sum(abs(bbeta))
  return(-ans)
}
