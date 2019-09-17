sub_daarem_lasso_gaussian_pgtn <- function(par, X, y, lambda, stplngth, nlag, a1, kappa, 
                                           maxiter, tol, mon.tol, cycl.mon.tol, sub.size, sub.type="random") {
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
    possible.index <- 1:num.params
    x.aa <- f.aa <- double(num.params)
    indf <- indx <- 1:num.params
    while(k < maxiter) {
        count <- count + 1

       # tmpf <- SoftThresh_ind(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) 
       # fnew <- tmpf$beta.hat - xnew
        fnew <- SoftThresh(xnew + stplngth*(Xty - crossprod(X,X%*%xnew)), lambda=lambda*stplngth) - xnew
        ss.resids <- sqrt(crossprod(fnew))

        
        if(ss.resids < tol & count==nlag) break

        np <- count
        if (count==1 & k > 20){
      #  if ( k%%10==0) {
             x.aa <- f.aa <- double(num.params)
            # ff.abs <- rowSums(abs(Fdiff))
            # xx.abs <- rowSums(abs(Xdiff))
            # ff.abs <- abs(fnew)
             ff.abs <- abs(fnew - fold) + abs(fnew)
             xx.abs <- abs(xnew - xold)
             indf <- which( ff.abs > .01*mean(ff.abs))
             indx <- which( xx.abs > .01*mean(xx.abs))
             #print(c(length(ind1), length(ind2), length(ind)))
        }
        Fdiff[,count] <- fnew - fold
        Xdiff[,count] <- xnew - xold
      
        Ftmp <- matrix(Fdiff[indf,1:np], nrow=length(indf), ncol=np)
        Xtmp <- matrix(Xdiff[indx,1:np], nrow=length(indx), ncol=np)
        fnew.tmp <- fnew[indf]
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
        
             x.aa[indx] <- Xtmp%*%gamma_vec
             f.aa[indf] <- Ftmp%*%gamma_vec

             xbar <- xnew - x.aa
             fbar <- fnew - f.aa

             x.propose <- xbar + fbar
             new.objective.val <- try(LassoObjFn_pgtn(x.propose, X, Xty, yty, lasso.pen), silent=TRUE)
     
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
       } else {
            ## Keep delta the same
            fold <- fnew
            xold <- xnew

            xnew <- fold + xold
      
            obj_funvals[k+2] <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)
     
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
    value.obj <- LassoObjFn_pgtn(xnew, X, Xty, yty, lasso.pen)

    if(k >= maxiter) {
        conv <- FALSE
        warning("Algorithm did not converge")
    }
    return(list(par=c(xnew), fpevals = k, value.objfn=value.obj, objfevals=obj.evals, 
                convergence=conv, objfn.track=obj_funvals, stplngth=stplngth))
}
