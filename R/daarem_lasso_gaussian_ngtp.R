daarem_lasso_gaussian_ngtp <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, 
                                       tol, mtol, cycl.mon.tol) {
    num.params <- ncol(X)
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
    obj_funvals <- rep(NA, maxiter + 2)

    XtX <- crossprod(X, X)
    yty <- sum(y*y)
    Xty <- crossprod(X, y)

    xold <- par
    xnew <- SoftThresh(xold + stplngth*(Xty - XtX%*%xold), lambda=lambda*stplngth)
    obj_funvals[1] <- LassoObjFn_ngtp(xold, XtX, Xty, yty, lasso.pen)
    obj_funvals[2] <- LassoObjFn_ngtp(xnew, XtX, Xty, yty, lasso.pen)

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
    if(length(mtol)==2) {
        mon.tol <- mtol[1]
    } else if(length(mtol)==1) {
        mon.tol <- mtol
    }
    while(k < maxiter) {
        count <- count + 1

        ## perform fixed-point update
        fnew <- SoftThresh(xnew + stplngth*(Xty - XtX%*%xnew), lambda=lambda*stplngth) - xnew
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

        Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
        tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
        lambda.ridge <- tmp_lam$lambda
        r.penalty <- tmp_lam$rr

        dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
        gamma_vec <- tmp$v%*%dd

        xbar <- xnew - drop(Xtmp%*%gamma_vec)
        fbar <- fnew - drop(Ftmp%*%gamma_vec)
        x.propose <- xbar + fbar
        new.objective.val <- try(LassoObjFn_ngtp(x.propose, XtX, Xty, yty, lasso.pen), silent=TRUE)
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

                obj_funvals[k+2] <- LassoObjFn_ngtp(xnew, XtX, Xty, yty, lasso.pen)
                obj.evals <- obj.evals + 1
                #num.em <- num.em + 1
            }
       } else {
           ## Keep delta the same
           fold <- fnew
           xold <- xnew
           xnew <- fold + xold

           obj_funvals[k+2] <- LassoObjFn_ngtp(xnew, XtX, Xty, yty, lasso.pen)  
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
            if(length(mtol)==2) {
                mon.tol= ifelse(mon.tol==mtol[1], mtol[2], mtol[1])
            } 
       }
       shrink.target <-  1/(1 + a1^(kappa - shrink.count))
       k <- k+1
    }
    obj_funvals <- obj_funvals[!is.na(obj_funvals)]
    value.obj <- LassoObjFn_ngtp(xnew, XtX, Xty, yty, lasso.pen)

    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals,
                convergence=conv, objfn.track=obj_funvals, stplngth=stplngth))
}
