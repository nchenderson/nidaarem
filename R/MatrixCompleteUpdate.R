MatrixCompleteUpdate <- function(par, lambda, Y, PY, ind) {
  tt <- 1
  B <- matrix(par, nrow=nrow(PY), ncol=ncol(PY))
  PB <- B
  PB[-ind] <- 0.0
  
  Atmp <- B + tt*(PY - PB)
  ans <- svd_soft_thresh(Atmp, lambda*tt)
  return(c(ans))
}