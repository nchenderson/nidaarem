sdaarem <- function(par, fixptfn, objfn, ..., control=list()) {

  control.default <- list(maxiter=2000, order=10, tol=1.e-08, mon.tol=0.01, cycl.mon.tol=0.0, kappa=25, alpha=1.2)
  namc <- names(control)
  if (!all(namc %in% names(control.default))) {
    stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
  }
  control <- modifyList(control.default, control)

  maxiter <- control$maxiter
  tol <- control$tol
  mon.tol <- control$mon.tol  ## monotonicity tolerance
  cycl.mon.tol <- control$cycl.mon.tol
  a1 <- control$alpha
  kappa <- control$kappa

  num.params <- length(par)
  nlag <- min(control$order, ceiling(num.params/2))

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
  #num.em <- 0  ## number of EM fallbacks
  ell.star <- obj_funvals[2]
  x.aa <- f.aa <- double(num.params)
  select <- 1:num.params
  while(k < maxiter) {
    count <- count + 1

    xtry <- fixptfn(xnew, ...)
    fnew <-  xtry - xnew
    ss.resids <- sqrt(crossprod(fnew))

   
    if(ss.resids < tol & count==nlag) break

    np <- count
    if (count==1 & k > 1){
         x.aa <- f.aa <- double(num.params)
         ff.abs <- rowSums(abs(Fdiff))
         xx.abs <- rowSums(abs(Xdiff))
         ind1 <- which( ff.abs > .1*mean(ff.abs))
         ind2 <- which( xx.abs > .1*mean(xx.abs))
         select <- union(ind1, ind2)
         #print(c(length(ind1), length(ind2), length(ind)))
    }
    print(c(length(select), nrow(Fdiff)))
    
    Fdiff[,count] <- fnew - fold
    Xdiff[,count] <- xnew - xold
    ff.abs <- rowSums(abs(Fdiff))
    xx.abs <- rowSums(abs(Xdiff))
  #  print(head(Fdiff))
  #  print(mean(xx.abs))
  #  print(head(Xdiff))
  #  print(mean(ff.abs))
        
    Ftmp <- matrix(Fdiff[select,1:np], nrow=length(select), ncol=np)
    Xtmp <- matrix(Xdiff[select,1:np], nrow=length(select), ncol=np)
    fnew.tmp <- fnew[select]
       # Ftmp.full <- matrix(Fdiff[,1:np], nrow=length(fnew), ncol=np)
       # Xtmp.full <- matrix(Xdiff[,1:np], nrow=length(fnew), ncol=np)

    tmp <- La.svd(Ftmp)
    dvec <- tmp$d
    uy <- crossprod(tmp$u, fnew.tmp)
    uy.sq <- uy*uy

    max.d <- max(tmp$d)
    min.d <- min(tmp$d)
    cond.number <- ifelse(max.d==min.d, 1, max.d/min.d)  ## to take care of cases where max.d=min.d=0
    if(cond.number > 1e10) {
            print('large condition number')
            shrink.count <- shrink.count - 2
    }

    ### Still need to compute Ftf
    Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew.tmp))^2))
    tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
    lambda.ridge <- tmp_lam$lambda
    r.penalty <- tmp_lam$rr

    dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
    gamma_vec <- crossprod(tmp$vt, dd)

    if(class(gamma_vec) != "try-error"){

      x.aa[select] <- Xtmp%*%gamma_vec
      f.aa[select] <- Ftmp%*%gamma_vec

      xbar <- xnew - x.aa
      fbar <- fnew - f.aa

      x.propose <- xbar + fbar
      new.objective.val <- try(objfn(x.propose, ...), silent=TRUE)
      obj.evals <- obj.evals + 1

      if(class(new.objective.val) != "try-error" & !is.na(obj_funvals[k+1]) &
         !is.nan(new.objective.val)) {
        if(new.objective.val >= obj_funvals[k+1] - mon.tol) {
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
          #num.em <- num.em + 1
        }
      } else {
        ## Keep delta the same
        fold <- fnew
        xold <- xnew

        xnew <- fold + xold
        obj_funvals[k+2] <- objfn(xnew, ...)
        obj.evals <- obj.evals + 1
        count <- 0
        #num.em <- num.em + 1
      }
    } else {
      ## Keep delta the same
      fold <- fnew
      xold <- xnew

      xnew <- fold + xold
      obj_funvals[k+2] <- objfn(xnew, ...)
      obj.evals <- obj.evals + 1
      count <- 0
      #num.em <- num.em + 1
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
    }
    
    shrink.target <-  1/(1 + a1^(kappa - shrink.count))
    k <- k+1
  }
  obj_funvals <- obj_funvals[!is.na(obj_funvals)]
  value.obj <- objfn(xnew, ...)
  if(k >= maxiter) {
    conv <- FALSE
    warning("Algorithm did not converge")
  }
  return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, convergence=conv, objfn.track=obj_funvals))
}

