nidaarem <- function(par, fixptfn, objfn, ..., control=list()) {
  
  if("objfn.check" %in% names(control)) {
    if(control$objfn.check) {
      control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=1, cycl.mon.tol=0.0, 
                              objfn.check=TRUE, kappa=25, alpha=1.2)
    }
    else if(!control$objfn.check) {
      control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=0.95, cycl.mon.tol=0.0, 
                              objfn.check=TRUE, kappa=25, alpha=1.2)
    }
  } else {
    control.default <- list(maxiter=2000, order=5, tol=1.e-08, mon.tol=1, cycl.mon.tol=0.0, 
                            objfn.check=TRUE, kappa=25, alpha=1.2)
  }
  control <- modifyList(control.default, control)
  
  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa
  resid.tol <- control$resid.tol
  
  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))
  
  if(!missing(objfn)) {
    par.init <- par
    neirun <- NesterovInitialize(par=par.init, fixptfn=fixptfn, objfn = objfn, 
                                 control=list(maxiter=maxiter, tol=tol), ...)
    
    nest.fpevals <- neirun$fpevals
    par <- neirun$par
    if(!neirun$convergence) {
        ans <- daaram_base_objfn(par, fixptfn, objfn, maxiter, tol, mon.tol, 
                                 cycl.mon.tol, a1, kappa, num.params, nlag, ...)
    } else {
        ans <- neirun
    }
  } else {
    ans <- daarem_base_noobjfn(par, fixptfn, maxiter, tol, resid.tol, 
                               a1, kappa, num.params, nlag, ...) 
  }
  if(!ans$convergence) {
    warning("Algorithm did not converge")
  }
  return(ans)
}