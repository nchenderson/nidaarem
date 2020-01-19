MatrixCompleteObj <- function(par, lambda, Y, PY, ind) {
  B <- matrix(par, nrow=nrow(PY), ncol=ncol(PY))
  svdB <- svd(B)
  ans <- sum((par[ind] - Y)*(par[ind] - Y))/2 + lambda*sum(svdB$d)
  return(-ans)
}