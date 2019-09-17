daarem_lasso_gaussian_pgtn2 <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, maxiter, tol, mon.tol, cycl.mon.tol) {
    num.params <- ncol(X)
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0.0, nrow=num.params, ncol=nlag)
    resid_vals <- objfn_track <- rep(NA, maxiter + 2)
    kr <- .01
    rho <- .95

    yty <- sum(y*y)
    Xty <- crossprod(X, y)
    
    xold <- par
    xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X, X%*%xold)), lambda=lambda*stplngth)
    fold <- xnew - xold
    objfn_track[1] <-  LassoObjFn_pgtn(xold, X, Xty, yty, lasso.pen)
    resid_vals[1] <- sqrt(crossprod(fold))
    fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
    resid_vals[2] <- sqrt(crossprod(fnew))
    ss.resids <- resid_vals[2]
    objfn_track[2] <-  LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
    krg <- ss.resids    
        
    k <- 1
    count <- 0
    shrink.count <- 0
    shrink.target <- 1/(1 + a1^kappa)
 
  
    lambda.ridge <- 100000
    r.penalty <- 0
    conv <- TRUE
    num.em <- 0  ## number of EM fallbacks
    n.aa <- 0
    DD <- 10
    n.viol <- 0
    while(k < maxiter) {
        count <- count + 1

        Fdiff[,count] <- fnew - fold
        Xdiff[,count] <- xnew - xold

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
        if(cond.number > 1e14) {
            shrink.count <- shrink.count - 2
        }

        ### Still need to compute Ftf
        Ftf <- sqrt(sum(as.vector(crossprod(Ftmp, fnew))^2))
        tmp_lam <- DampingFind(uy.sq, dvec, a1, kappa, shrink.count, Ftf, lambda.start=lambda.ridge, r.start=r.penalty)
        lambda.ridge <- tmp_lam$lambda
        r.penalty <- tmp_lam$rr

        dd <- (dvec*uy)/(dvec^2 + lambda.ridge)
        gamma_vec <- crossprod(tmp$vt, dd)

        if(class(gamma_vec) != "try-error"){

             xbar <- xnew - drop(Xtmp%*%gamma_vec)
             fbar <- fnew - drop(Ftmp%*%gamma_vec)

             x.propose <- xbar + fbar
           
             ftmp <- SoftThresh(x.propose + stplngth*(Xty - crossprod(X,X%*%x.propose)), lambda=lambda*stplngth) - x.propose
             ss.tmp <- sqrt(crossprod(ftmp))
             
             theor.val <- (DD*krg)/((n.aa + 1)^(1 + .05))
             #print(c(ss.tmp, theor.val))
            # if(ss.tmp > theor.val) {
            #      n.viol <- n.viol + 1
            # }
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

                    xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X,X%*%xold)), lambda=lambda*stplngth)
                    fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
                    ss.resids <- sqrt(crossprod(fnew))
             }
             
       } else {
            ## Keep delta the same
            fold <- fnew
            xold <- xnew

            xnew <- SoftThresh(xold + stplngth*(Xty - crossprod(X,X%*%xold)), lambda=lambda*stplngth)
            fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
            ss.resids <- sqrt(crossprod(fnew))
            
            count <- 0
       }
       if(ss.resids < tol & count==nlag) break
      
       resid_vals[k + 2] <- ss.resids
       objfn_track[k + 2] <-  LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
       if(count==nlag) {
            count <- 0
            ## restart count
       }
       
       shrink.target <-  1/(1 + a1^(kappa - shrink.count))
       k <- k+1
    }
    value.obj <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
    resid_vals <- resid_vals[!is.na(resid_vals)]
    objfn_track <- objfn_track[!is.na(objfn_track)]
    
    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfn.track=objfn_track,
                convergence=conv, resid.track=resid_vals, stplngth=stplngth, num.aa=n.aa, n.viol=n.viol))
}
