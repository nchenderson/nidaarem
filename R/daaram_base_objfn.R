daaram_base_objfn <- function(par, fixptfn, objfn, maxiter, tol, mtol, 
                              cycl.mon.tol, a1, kappa, num.params, nlag, ...) {
  
  Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
  obj_funvals <- rep(NA, maxiter + 2)
  
  xold <- par
  xnew <- fixptfn(xold, ...)
  obj_funvals[1] <- objfn(xold, ...)
  obj_funvals[2] <- objfn(xnew, ...)
  likchg <- obj_funvals[2] - obj_funvals[1]
  obj.evals <- 2
  
  fold <- xnew - xold
  k <- 1
  count <- 0
  shrink.count <- 0
  shrink.target <- 1/(1 + a1^kappa)
  lambda.ridge <- 100000
  r.penalty <- 0
  conv <- TRUE
  ell.star <- obj_funvals[2]
  if(length(mtol)==2) {
    mon.tol <- mtol[1]
  } else if(length(mtol)==1) {
    mon.tol <- mtol
  }
  while(k < maxiter) {
    count <- count + 1
    
    fnew <- fixptfn(xnew, ...) - xnew
    ss.resids <- sqrt(crossprod(fnew))
    if(ss.resids < tol & count==nlag) break
    
    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    
    np <- count
    if(np==1) {
      Ftmp <- matrix(Fdiff[,1], nrow=num.params, ncol=np)
      Xtmp <- matrix(Xdiff[,1], nrow=num.params, ncol=np)  ## is this matrix function needed?
      
      #pow <- kappa - shrink.count
      #s.delta <- exp(-0.5*log1p(a1^pow))
      #dsq <- sum(Fdiff[,count]*Fdiff[,count])
      #dsq <- sum(Ftmp^2)
      #lambda.ridge <- dsq*(1/s.delta - 1)
      ## What if dsq = 0?  
      #gamma_vec <- sum(Fdiff[,count]*fnew)/(dsq + lambda.ridge)  
      #r.penalty = 0
    } else {
      Ftmp <- Fdiff[,1:np]
      Xtmp <- Xdiff[,1:np]  
    }
    tmp <- La.svd(Ftmp)
    dvec <- tmp$d
    dvec.sq <- dvec*dvec
    uy <- crossprod(tmp$u, fnew)
    uy.sq <- uy*uy
    
    ### Still need to compute Ftf
    Ftf <- sqrt(sum(uy.sq*dvec.sq))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr
    dd <- (dvec*uy)/(dvec.sq + lambda.ridge)
    gamma_vec <- crossprod(tmp$vt, dd)
    
    xbar <- xnew - drop(Xtmp%*%gamma_vec)
    fbar <- fnew - drop(Ftmp%*%gamma_vec)
    
    x.propose <- xbar + fbar
    new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)
    obj.evals <- obj.evals + 1
    
    if(class(new.objective.val)[1] != "try-error" & !is.na(obj_funvals[k+1]) &
       !is.nan(new.objective.val)) {
      if(new.objective.val >= obj_funvals[k+1] - mon.tol) {  ## just change this line in daarem_base_noobjfn
        ## Increase delta
        obj_funvals[k+2] <- new.objective.val
        fold <- fnew
        xold <- xnew
        
        xnew <- x.propose
        shrink.count <- shrink.count + 1
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew
        
        xnew <- fold + xold
        obj_funvals[k+2] <- objfn(xnew, ...)
        obj.evals <- obj.evals + 1
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew
      
      xnew <- fold + xold
      obj_funvals[k+2] <- objfn(xnew, ...)
      obj.evals <- obj.evals + 1
      count <- 0
    }
    if(count==nlag) {
      count <- 0
      ## restart count
      ## make comparison here l.star vs. obj_funvals[k+2]
      if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
        ## Decrease delta
        shrink.count <- max(shrink.count - nlag, -2*kappa)
      }
      ell.star <- obj_funvals[k+2]
      if(length(mtol)==2) {
        mon.tol= ifelse(mon.tol==mtol[1], mtol[2], mtol[1])
      } 
    }
    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k+1
  }
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  value.obj <- objfn(xnew, ...)
  if(k >= maxiter) {
    conv <- FALSE
  }
  return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, convergence=conv, objfn.track=obj_funvals))
}
