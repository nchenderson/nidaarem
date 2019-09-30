
Nesterov <- function(par, fixptfn, objfn=NULL, restart=FALSE, ...,
                     control=list(maxiter=1000, tol=1e-7)) {

  maxiter <- control$maxiter
  tol <- control$tol
  ## 4 possible combinations of restart and objective function.
  if(is.null(objfn)) {
    base_fn <- "nest_no_obj"
  } else if(!is.null(objfn) & !restart) {
    base_fn <- "nest_obj_no_restart"
  } else if(!is.null(objfn) & restart) {
    base_fn <- "nest_obj_restart"
  } else if(is.null(objfn) & restart) {
    base_fn <- "nest_noobj_restart"
  }

  ans <- switch(base_fn, nest_no_obj = NesterovNoObjNR(par, fixptfn, ..., maxiter=maxiter, tol=tol),
                nest_obj_no_restart = NesterovObjNR(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol),
                nest_obj_restart = NesterovObjR(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol),
                nest_noobj_restart = NesterovObjNR(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol))
  return(ans)
}


NesterovObjNR <- function(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol) {
  pp <- length(par)
  conv <- TRUE
  beta.vecy <- beta.old <- par
  objfn.val <- rep(NA, maxiter+1)
  alpha <- 1
  objfn.val[1] <- objfn(beta.old, ...)
  for(k in 1:maxiter) {
    beta.new <- fixptfn(beta.vecy, ...)
    ## Update tseq
    alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
    mom <- (alpha - 1)/alpha.new
    alpha <- alpha.new
    beta.vecy <- beta.new + mom*(beta.new - beta.old)

    ss.resids <- sqrt(crossprod(beta.new - beta.old))
    beta.old <- beta.new
    objfn.val[k+1] <- objfn(beta.new, ...)
    if(ss.resids < tol) break
  }
  objfn.val <- objfn.val[!is.na(objfn.val)]
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  value.obj <- objfn.val[k+1]
  return(list(par=beta.old, value.objfn=value.obj, objfn.track=objfn.val, fpevals=k, convergence=conv))
}


NesterovObjR <- function(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol) {
  pp <- length(par)
  conv <- TRUE
  beta.vecy <- beta.old <- par
  objfn.val <- rep(NA, maxiter+1)
  alpha <- 1
  objfn.val[1] <- objfn(beta.old, ...)
  for(k in 1:maxiter) {
    beta.new <- fixptfn(beta.vecy, ...)
    oval <- objfn(beta.new, ...)
    if(oval < objfn.val[k]) {
      ## restart
      alpha <- 1
      beta.vecy <- beta.new
    } else {
      ## Update tseq
      alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
      mom <- (alpha - 1)/alpha.new
      alpha <- alpha.new
      beta.vecy <- beta.new + mom*(beta.new - beta.old)

      ss.resids <- sqrt(crossprod(beta.new - beta.old))
      if(ss.resids < tol) break
    }
    objfn.val[k+1] <- objfn(beta.new, ...)
    beta.old <- beta.new
  }
  objfn.val <- objfn.val[!is.na(objfn.val)]
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  value.obj <- objfn.val[k]
  return(list(par=beta.new, value.objfn=value.obj, objfn.track=objfn.val, fpevals=k, convergence=conv))
}


NesterovNoObjNR <- function(par, fixptfn, ..., maxiter=maxiter, tol=tol) {
  pp <- length(par)
  conv <- TRUE
  beta.vecy <- beta.old <- par
  alpha <- 1
  for(k in 1:maxiter) {
    beta.new <- fixptfn(beta.vecy, ...)
    ## Update tseq
    alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
    mom <- (alpha - 1)/alpha.new
    alpha <- alpha.new
    beta.vecy <- beta.new + mom*(beta.new - beta.old)

    ss.resids <- sqrt(crossprod(beta.new - beta.old))
    beta.old <- beta.new
    if(ss.resids < tol) break
  }
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  return(list(par=beta.old, fpevals=k, convergence=conv))
}


NesterovNoObjR <- function(par, fixptfn, objfn, ..., maxiter=maxiter, tol=tol) {
  pp <- length(par)
  conv <- TRUE
  beta.vecy <- beta.old <- par
  objfn.val <- rep(NA, maxiter+1)
  alpha <- 1
  objfn.val[1] <- objfn(beta.old, ...)
  for(k in 1:maxiter) {
    beta.new <- fixptfn(beta.vecy, ...)
    gc <- sum((beta.vecy - beta.new)*(beta.new - beta.old))
    if(gc > 0) {
      ## restart
      alpha <- 1
      beta.vecy <- beta.new
    } else {
      ## Update tseq
      alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
      mom <- (alpha - 1)/alpha.new
      alpha <- alpha.new
      beta.vecy <- beta.new + mom*(beta.new - beta.old)

      ss.resids <- sqrt(crossprod(beta.new - beta.old))
      if(ss.resids < tol) break
    }
    objfn.val[k+1] <- objfn(beta.new, ...)
    beta.old <- beta.new
  }
  objfn.val <- objfn.val[!is.na(objfn.val)]
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  value.obj <- objfn.val[k]
  return(list(par=beta.new, fpevals=k, conv=convergence))
}
