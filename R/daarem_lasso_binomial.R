daarem_lasso_binomial <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol) {
    num.params <- ncol(X)
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
    obj_funvals <- rep(NA, maxiter + 2)

    Xty <- crossprod(X, y)
    
    xold <- par
    phat <- expit(X%*%xold)
    xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X, phat)), lambda=lambda*stplngth)
    obj_funvals[1] <- LogisticObjFn(xold, X, Xty, lasso.pen)
    obj_funvals[2] <- LogisticObjFn(xnew, X, Xty, lasso.pen)
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
    num.em <- 0  ## number of EM fallbacks
    ell.star <- obj_funvals[2]
    while(k < maxiter) {
        count <- count + 1

        rr <- y - expit(X%*%xnew) ## residuals
        fnew <- SoftThresh(xnew + stplngth*crossprod(X, rr), lambda=lambda*stplngth) - xnew
        ss.resids <- sqrt(crossprod(fnew))

        Fdiff[,count] <- fnew - fold
        Xdiff[,count] <- xnew - xold
        if(ss.resids < tol & count==nlag) break


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

        if(class(gamma_vec) != "try-error"){

             xbar <- xnew - drop(Xtmp%*%gamma_vec)
             fbar <- fnew - drop(Ftmp%*%gamma_vec)

             x.propose <- xbar + fbar
             new.objective.val <- try(LogisticObjFn(x.propose, X, Xty, lasso.pen), silent=TRUE)
     
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
                    
                    ### Do we need to re-compute everything if we fall back?
                    obj_funvals[k+2] <- LogisticObjFn(xnew, X, Xty, lasso.pen)
                    obj.evals <- obj.evals + 1
                    #num.em <- num.em + 1
                 }
             } else {
                 ## Keep delta the same
                 fold <- fnew
                 xold <- xnew

                 xnew <- fold + xold

                 obj_funvals[k+2] <- LogisticObjFn(xnew, X, Xty, lasso.pen)  ### need to add ngtp here?
       
                 obj.evals <- obj.evals + 1
                 count <- 0
                 #num.em <- num.em + 1
            }
       } else {
            ## Keep delta the same
            fold <- fnew
            xold <- xnew

            xnew <- fold + xold
      
            obj_funvals[k+2] <- LogisticObjFn(xnew, X, Xty, lasso.pen)
     
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
    value.obj <- LogisticObjFn(xnew, X, Xty, lasso.pen)

    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, 
                convergence=conv, objfn.track=obj_funvals, stplngth=stplngth))
}

