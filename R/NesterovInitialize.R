NesterovInitialize <-  function(par, fixptfn, objfn=NULL, test="monotone", ...,
                                control=list(maxiter=500, tol=1e-7)) {
    control.default <- list(maxiter=500, tol=1e-7)
    namc <- names(control)
    if (!all(namc %in% names(control.default))) {
       stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
    }
    control <- modifyList(control.default, control)
    maxiter <- control$maxiter
    tols <- control$tol
    
    pp <- length(par)
    conv <- TRUE
    beta.vecy <- beta.old <- par
    objfn.val <- rep(NA, maxiter+1)
    alpha <- 1
    objfn.val[1] <- objfn(beta.old, ...)
    conv <- FALSE
    for(k in 1:maxiter) {
        beta.new <- fixptfn(beta.vecy, ...) ## x_k, y_k = beta.vecy
        oval <- objfn(beta.new, ...)
        if(test=="monotone") {
           done1 <- oval < objfn.val[k] 
        } else if(test == "gradient") {
           done1 <- sum((beta.vecy - beta.new)*(beta.new - beta.old)) > 0 
        }
        ### perform a second "done test" here.
        #if(done1) {
        #    xold <- beta.new  ##x0
        #    xnew <- fixptfn(xold, ...) ## x1
        #    fold <- xnew - xold ## f0 and Delta x0
        #    fnew <- fixptfn(xnew, ...) - xnew ## f1
        #    deltaf <- fnew - fold ## Delta f0
        #    gamma.aa <- sum(deltaf*f1)/sum(deltaf*deltaf)
        #    x.aa <- xnew + fnew - gamma.aa*fnew
        #    oval.aa <- objfn(x.aa, ...)
        #    aa.test <- ifelse(oval.aa > , TRUE, FALSE)
        #}
        if(done1) {
            break
        } 
       
        ## Update tseq
        alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
        mom <- (alpha - 1)/alpha.new
        alpha <- alpha.new
        
        beta.vecy <- beta.new + mom*(beta.new - beta.old)
        objfn.val[k+1] <- objfn(beta.new, ...)
        
        ss.resids <- sqrt(crossprod(beta.new - beta.old))
        done2 <- ss.resids < tols
        if(done2) {
           conv <- TRUE
           break
        }
        beta.old <- beta.new
    }
    objfn.val <- objfn.val[!is.na(objfn.val)]
    value.obj <- objfn.val[k]
    return(list(par=beta.new, value.objfn=value.obj, objfn.track=objfn.val, fpevals=k, convergence=conv))
}
