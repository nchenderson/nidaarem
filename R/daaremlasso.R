daarem.lasso <- function(par, X, y, lambda, stplngth=NULL, nesterov.init=FALSE,
                         family = c("gaussian", "binomial"), control=list()) {

  if("objfn.check" %in% names(control)) {
      if(control$objfn.check) {
           control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=1, cycl.mon.tol=0.0, 
                                  objfn.check=TRUE, kappa=25, alpha=1.2)
      }
      else if(!control$objfn.check) {
           control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=0.95, cycl.mon.tol=0.0, 
                                  objfn.check=TRUE, kappa=25, alpha=1.2)
      }
  } else {
      control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=1, cycl.mon.tol=0.0, 
                              objfn.check=TRUE, kappa=25, alpha=1.2)
  }
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

  if(family=="gaussian" & control$objfn.check) {
      base_fn <- ifelse(n > p, "gauss_ngtp", "gauss_pgtn")
  } else if(family=="gaussian" & !control$objfn.check) {
      base_fn <- ifelse(n > p, "gauss_ngtp2", "gauss_pgtn2")
  } else if(family=="binomial" & control$objfn.check) {
      base_fn <- "binomial_b"
  } else if(family == "binomial" & !control$objfn.check) {
      base_fn <- "binomial_b2"
  }
  #if(nesterov.init=TRUE) {
  #    tmp <- NesterovInitialize(par, fixptfn, )
  #    par <- tmp$par
  #}
  nest.fpevals <- 0
  if(nesterov.init & family=="gaussian") {
      par.init <- par
      neirun <- NesterovInitialize(par=par.init, fixptfn=GDLassoStep, objfn = LassoObjFn, 
                                   test = "monotone", X=X, y=y, lambda=lambda, stplngth=stplngth, 
                                   control=list(maxiter=maxiter, tol=tol))
      nest.fpevals <- neirun$fpevals
      par <- neirun$par
      if(!neirun$convergence) {
          ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_pgtn = daarem_lasso_gaussian_pgtn(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        binomial_b = daarem_lasso_binomial(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        binomial_b2 = daarem_lasso_binomial2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol))
          ans$fpevals <- ans$fpevals + nest.fpevals
          ans$objfn.track <- c(neirun$objfn.track, ans$objfn.track)
      } else {
          ans <- neirun
      }
  } else if(nesterov.init & family=="binomial") {
      Xty <- crossprod(X, y)
      par.init <- par
      neirun <- NesterovInitialize(par=par.init, fixptfn=GDLogisticStep, objfn = LogisticObjFn, 
                                 test = "monotone", X=X, Xty=Xty, lambda=lambda, stplngth=stplngth, 
                                 control=list(maxiter=maxiter, tol=tol))
      nest.fpevals <- neirun$fpevals
      par <- neirun$par
      if(!neirun$convergence) {
          ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_pgtn = daarem_lasso_gaussian_pgtn(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        binomial_b = daarem_lasso_binomial(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
                        binomial_b2 = daarem_lasso_binomial2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol))
          ans$fpevals <- ans$fpevals + nest.fpevals
           ans$objfn.track <- c(neirun$objfn.track, ans$objfn.track)
      } else {
          ans <- neirun
      }
  } else {
     ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
         gauss_pgtn = daarem_lasso_gaussian_pgtn(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
         gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
         gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
         binomial_b = daarem_lasso_binomial(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol),
         binomial_b2 = daarem_lasso_binomial2(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter - nest.fpevals, tol, mon.tol, cycl.mon.tol))
  }
  #ans$fpevals <- ans$fpevals + nest.fpevals
  #ans$nest.iters <- nest.fpevals
  return(ans)
}

