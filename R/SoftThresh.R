SoftThresh <- function(x, lambda) {
  x_gelambda <- which(x >= lambda)
  x_lelambda <- which(x <= -lambda)
  sol <- double(length(x))
  sol[x_gelambda] <- x[x_gelambda] - lambda
  sol[x_lelambda] <- x[x_lelambda] + lambda
  return(sol)
}