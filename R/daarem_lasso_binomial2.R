daarem_lasso_binomial2 <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, 
                                   tol, mtol, cycl.mon.tol) {
    num.params <- ncol(X)
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
    resid_vals <- objfn_track <- rep(NA, maxiter + 2)
    kr <- .01
    if(length(mtol)==2) {
      rho <- .95
    } else {
      rho <- mtol[1]
    }
    
    Xty <- crossprod(X, y)
    
    xold <- par
    phat <- expit(X%*%xold)
    xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X, phat)), lambda=lambda*stplngth) 
    fold <- xnew - xold
    objfn_track[1] <-  LogisticObjFn(xold, X, Xty, lasso.pen)
    resid_vals[1] <- sqrt(crossprod(fold))
    phat <- expit(X%*%xnew)
    fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X, phat)), lambda=lambda*stplngth) - xnew
    resid_vals[2] <- sqrt(crossprod(fnew))
    ss.resids <- resid_vals[2]
    objfn_track[2] <-  LogisticObjFn(xnew, X, Xty, lasso.pen)
    krg <- ss.resids 

    k <- 1
    count <- 0
    shrink.count <- 0
    shrink.target <- 1/(1 + a1^kappa)
  
    lambda.ridge <- 100000
    r.penalty <- 0
    conv <- TRUE
    fallback = 0
    n.aa <- 0
    DD <- 10
    while(k < maxiter) {
        count <- count + 1

        #rr <- y - expit(X%*%xnew) ## residuals
        #fnew <- SoftThresh(xnew + stplngth*crossprod(X, rr), lambda=lambda*stplngth) - xnew
        #ss.resids <- sqrt(crossprod(fnew))

        Fdiff[,count] <- fnew - fold
        Xdiff[,count] <- xnew - xold

        np <- count
        Ftmp <- matrix(Fdiff[,1:np], nrow=length(fnew), ncol=np)
        Xtmp <- matrix(Xdiff[,1:np], nrow=length(fnew), ncol=np)

        tmp <- svd(Ftmp)
        dvec <- tmp$d
        uy <- crossprod(tmp$u, fnew)
        uy.sq <- uy*uy

        max.d <- max(tmp$d)
        min.d <- min(tmp$d)
        cond.number <- ifelse(max.d==min.d, 1, max.d/min.d)  ## to take care of cases where max.d=min.d=0
        if(cond.number > 1e14) {
            shrink.count <- shrink.count - 2
        }

        ### Still need to compute Ftf
        Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
        tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
        lambda.ridge <- tmp_lam$lambda
        r.penalty <- tmp_lam$rr

        dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
        gamma_vec <- tmp$v%*%dd

        xbar <- xnew - drop(Xtmp%*%gamma_vec)
        fbar <- fnew - drop(Ftmp%*%gamma_vec)
        x.propose <- xbar + fbar
           
        rr <- y - expit(X%*%x.propose)
        ftmp <- SoftThresh(x.propose + stplngth*crossprod(X, rr), lambda=lambda*stplngth) - x.propose
        ss.tmp <- sqrt(crossprod(ftmp))
             
        theor.val <- (DD*krg)/((n.aa + 1)^(1 + .05))
        if(ss.tmp <= min(ss.resids + rho^k, theor.val)) {
                 ## Increase delta
                    fold <- fnew
                    xold <- xnew
                     
                    xnew <- x.propose
                    fnew <- ftmp
                    ss.resids <- ss.tmp
                    
                    shrink.count <- shrink.count + 1
                    n.aa <- n.aa + 1
         } else {
                 ## Keep delta the same
                    fold <- fnew
                    xold <- xnew
                    
                    rr.old <- y - expit(X%*%xold)
                    xnew <- SoftThresh(xold + stplngth*crossprod(X, rr.old), lambda=lambda*stplngth)
                    rr.new <- y - expit(X%*%xnew)
                    fnew <- SoftThresh(xnew + stplngth*crossprod(X, rr.new), lambda=lambda*stplngth) - xnew
                    ss.resids <- sqrt(crossprod(fnew))
                    fallback <- fallback + 1
          }
          resid_vals[k + 2] <- ss.resids
          objfn_track[k + 2] <-  LogisticObjFn(xnew, X, Xty, lasso.pen)
          if(ss.resids < tol) break
          if(count==nlag) {
            count <- 0
            if(length(mtol)==2) {
              rho = ifelse(rho==mtol[1], mtol[2], mtol[1])
            } 
          }
          shrink.target <-  1/(1 + a1^(kappa - shrink.count))
          k <- k+1
    }
    objfn_track <- objfn_track[!is.na(objfn_track)]
    value.obj <- LogisticObjFn(xnew, X, Xty, lasso.pen)

    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfn.track=objfn_track,
                convergence=conv, resid.track=resid_vals, stplngth=stplngth, fallback=fallback))
}

