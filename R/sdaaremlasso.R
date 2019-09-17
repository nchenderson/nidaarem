sdaarem.lasso <- function(par, X, y, lambda, stplngth=NULL, 
                         family = c("gaussian", "binomial"), sub.size=500, sub.type="random", control=list()) {

  control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, alpha=1.2)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)

  family = match.arg(family)

  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotone tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa

  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

  lasso.pen <- lambda
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)

  ### Compute step length
  if(is.null(stplngth) & family=="gaussian") {
      stplngth <- FindStepLengthGaussian(X, y, lambda)
      ## Migth be better to modify FindStepLengthGaussian so that
      ## it takes XtX as an argument?
  } else if(is.null(stplngth) & family=="binomial") {
      stplngth <- FindStepLengthBinomial(X, y, lambda)
  }

  if(family=="gaussian" & sub.type=="random") {
      base_fn <- ifelse(n > p, "gauss_ngtp", "gauss_pgtn")
  } else if(family == "binomial") {
      base_fn <- "binomial_b"
  } else if(family=="gaussian" & sub.type=="bb") {
      base_fn <- "gauss_pgtn2"
  }
  ans <- switch(base_fn, gauss_ngtp = sub_daarem_lasso_gaussian_ngtp(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol, sub.size, sub.type),
         gauss_pgtn = sub_daarem_lasso_gaussian_pgtn(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol, sub.size, sub.type),
         binomial_b = sub_daarem_lasso_binomial(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol, sub.size),
         gauss_pgtn2 = sub_daarem_lasso_gaussian_pgtn(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol, sub.size, sub.type))
  return(ans)
}



