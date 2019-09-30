GDLassoStep <- function(par, X, y, lambda, stplngth) {
  beta.vec <- par
  X.beta <- X%*%beta.vec
  beta.vec <- SoftThresh(beta.vec + stplngth*crossprod(X, y - X.beta), 
                         lambda=lambda*stplngth)
  return(beta.vec)
}