svd_soft_thresh <- function(A, lambda) {
  tmp <- svd(A)
  d_thresh <- pmax(tmp$d - lambda, 0)
  ans <- tcrossprod(t(t(tmp$u)*d_thresh),tmp$v)
  return(ans)
}