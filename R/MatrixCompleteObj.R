MatrixCompleteObj <- function(par, lambda, A, PA, ind) {
  B <- matrix(par, nrow=nrow(PA), ncol=ncol(PA))
  svdB <- svd(B)
  
  ans <- sum((par[ind] - A)*(par[ind] - A))/2 + lambda*sum(svdB$d)
  return(-ans)
}