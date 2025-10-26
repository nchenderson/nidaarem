MatrixCompleteUpdate <- function(par, lambda, A, PA, ind) {
  ## ind is the vector indicating the observed indices
  tt <- 1
  B <- matrix(par, nrow=nrow(PA), ncol=ncol(PA))
  PB <- B
  PB[-ind] <- 0.0
  
  Atmp <- B + tt*(PA - PB)
  ans <- svd_soft_thresh(Atmp, lambda*tt)
  return(c(ans))
}

