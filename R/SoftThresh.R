SoftThresh <- function(x, lambda) {
  sign(x)*pmax(abs(x) - lambda, 0.0)
}
