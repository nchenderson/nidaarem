daarem_lasso_gaussian_pgtn <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, 
                                       tol, mtol, cycl.mon.tol) {
    num.params <- ncol(X)
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
    obj_funvals <- rep(NA, maxiter + 2)

    yty <- sum(y*y)
    Xty <- crossprod(X, y)
    
    xold <- par
    xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X, X%*%xold)), lambda=lambda*stplngth)
    obj_funvals[1] <- LassoObjFn_pgtn(xold, X, Xty, yty, lasso.pen)
    obj_funvals[2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
    ftmp <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
    ss.resids <- sqrt(crossprod(ftmp))
    
    likchg <- obj_funvals[2] - obj_funvals[1]
    obj.evals <- 2

    fold <- xnew - xold
    k <- 1
    count <- 0
    shrink.count <- 0
    shrink.target <- 1/(1 + a1^kappa)
    n.aa <- 0
    DD <- 100
    krg <- ss.resids    
    alpha <- 1
  
    lambda.ridge <- 100000
    r.penalty <- 0
    conv <- TRUE
    num.em <- 0  ## number of EM fallbacks
    ell.star <- obj_funvals[2]
    pred.val <- actual.val <- rep(NA, ceiling(maxiter + 2) + 5)
    h <- 1
    if(length(mtol)==2) {
        mon.tol <- mtol[1]
    } else if(length(mtol)==1) {
        mon.tol <- mtol
    }
    while(k < maxiter) {
        count <- count + 1

        fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
        ss.resids <- sqrt(crossprod(fnew))

        Fdiff[,count] <- fnew - fold
        Xdiff[,count] <- xnew - xold
        if(ss.resids < tol & count==nlag) break


        np <- count
        Ftmp <- matrix(Fdiff[,1:np], nrow=length(fnew), ncol=np)
        Xtmp <- matrix(Xdiff[,1:np], nrow=length(fnew), ncol=np)

        tmp <- La.svd(Ftmp)
        dvec <- tmp$d
        uy <- crossprod(tmp$u, fnew)
        uy.sq <- uy*uy

        max.d <- max(tmp$d)
        min.d <- min(tmp$d)
        cond.number <- ifelse(max.d==min.d, 1, max.d/min.d)  ## to take care of cases where max.d=min.d=0
        rho.max <- 1e16  ## try 1e16 instead?
        if(cond.number > sqrt(rho.max)) {
            shrink.count <- shrink.count - 1
        }

        ### Still need to compute Ftf
        Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
        tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
        lambda.ridge <- tmp_lam$lambda
        r.penalty <- tmp_lam$rr
        
        ##################################################
        lambda.tmp <- -1e-16
        if(min.d*min.d*rho.max - max.d*max.d < 0) { ## if "ridge regression condition number (with lambda=0) is greater than rho.max"
           lambda.tmp <- (rho.max*min.d*min.d - max.d*max.d)/(1 - rho.max)
        }
        #if(lambda.tmp > lambda.ridge) {
        #    print(c(max.d, min.d, min.d*min.d*rho.max, max.d*max.d, lambda.ridge))
        #}
        lambda.ridge <- max(c(lambda.tmp, lambda.ridge, 0.0))
        ####################################################  
        
        dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
        gamma_vec <- crossprod(tmp$vt, dd)

        xbar <- xnew - drop(Xtmp%*%gamma_vec)
        fbar <- fnew - drop(Ftmp%*%gamma_vec)
        x.propose <- xbar + fbar
        new.objective.val <- try(LassoObjFn_pgtn(x.propose, X, Xty, yty, lasso.pen), silent=TRUE)
        obj.evals <- obj.evals + 1

        ftmp <- SoftThresh(x.propose + stplngth*(Xty - crossprod(X,X%*%x.propose)), lambda=lambda*stplngth) - x.propose
        ss.tmp <- sqrt(crossprod(ftmp))
             
        theor.val <- (DD*krg)/((n.aa + 1)^(1 + .05))
             
        if(class(new.objective.val) != "try-error" & !is.na(obj_funvals[k+1]) &
            !is.nan(new.objective.val)) {
             if(new.objective.val >= obj_funvals[k+1] - mon.tol & ss.tmp <= theor.val) {
                 ## Increase delta
                 obj_funvals[k+2] <- new.objective.val
                 fold <- fnew
                 xold <- xnew

                 xnew <- x.propose
                 shrink.count <- shrink.count + 1
                 n.aa <- n.aa + 1
             } else {
                 ## Keep delta the same
                 fold <- fnew
                 xold <- xnew
                 xnew <- fold + xold
                    
                 ## Update tseq
                 #if(new.objective.val < obj_funvals[k+1]) {
                 #    alpha <- 1
                 #}
                 #alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
                 #mom <- (alpha - 1)/alpha.new
                 #alpha <- alpha.new
                 #xnew <- xnew + mom*(xnew - xold)
                    
                 ### Do we need to re-compute everything if we fall back?
                 obj_funvals[k+2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
                 obj.evals <- obj.evals + 1
                 #num.em <- num.em + 1
             }
        } else {
             ## Keep delta the same
             fold <- fnew
             xold <- xnew
             xnew <- fold + xold

             obj_funvals[k+2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)  ### need to add ngtp here?
             obj.evals <- obj.evals + 1
             count <- 0
             #num.em <- num.em + 1
       }
       if(count==nlag) {
            count <- 0
            ## restart count
            ## make comparison here l.star vs. obj_funvals[k+2]
            ## switch monotonicity parameter here.
            if(obj_funvals[k+2] < ell.star - cycl.mon.tol) {
                ## Decrease delta
               shrink.count <- max(shrink.count - nlag, -2*kappa)
            }
            if(length(mtol)==2) {
                mon.tol= ifelse(mon.tol==mtol[1], mtol[2], mtol[1])
            } 
            ell.star <- obj_funvals[k+2]
            xtmp <- SoftThresh(xnew + stplngth*(Xty - crossprod(X, X%*%xnew)), lambda=lambda*stplngth)
            objtmp <- LassoObjFn_pgtn(xtmp, X, Xty, yty, lasso.pen)
            pred.val[h] <- 2*nlag*(objtmp - ell.star) + ell.star
            actual.val[h] <- ell.star
              
            h <- h+1
            #################################################################
            #xold <- xnew
            #xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X, X%*%xold)), lambda=lambda*stplngth)
            ## need to add objective function value here
            #k <- k + 1
            #obj_funvals[k+2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
            #fold <- xnew - xold
            #################################################################
            ## xnew is x0 and y1
            #xold <- xnew
            #beta.new <- SoftThresh(xold + stplngth*(Xty - crossprod(X,X%*%xold)), lambda=lambda*stplngth) ##x1
            #fold <- beta.new - xold
            #k <- k+1
            #alpha <- 1
            #alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
            #mom <- (alpha - 1)/alpha.new
            #alpha <- alpha.new
            #xnew <- beta.new + mom*(beta.new - xold) ##y2
            #obj_funvals[k+2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
            
       }
       shrink.target <-  1/(1 + a1^(kappa - shrink.count))
       k <- k+1
    }
    obj_funvals <- obj_funvals[!is.na(obj_funvals)]
    value.obj <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)

    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, 
                convergence=conv, objfn.track=obj_funvals, stplngth=stplngth, n.aa=n.aa, pred.val=pred.val,
                actual.val=actual.val))
}
