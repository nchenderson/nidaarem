NesterovInitialize <-  function(par, fixptfn, objfn=NULL, test="monotone", ...,
                                control=list(maxiter=500)) {
    control.default <- list(maxiter=500)
    namc <- names(control)
    if (!all(namc %in% names(control.default))) {
       stop("unknown names in control: ", namc[!(namc %in% names(control.default))])
    }
    control <- modifyList(control.default, control)
    maxiter <- control$maxiter
    
    
    pp <- length(par)
    conv <- TRUE
    beta.vecy <- beta.old <- par
    objfn.val <- rep(NA, maxiter+1)
    alpha <- 1
    objfn.val[1] <- objfn(beta.old, ...)
    for(k in 1:maxiter) {
        beta.new <- fixptfn(beta.vecy, ...) ## x_k, y_k = beta.vecy
        oval <- objfn(beta.new, ...)
        if(test=="monotone") {
           done <- oval < objfn.val[k] 
        } else if(test == "gradient") {
           done <- sum((beta.vecy - beta.new)*(beta.new - beta.old)) > 0 
        }        
        if(done) {
            break
        } 
        ## Update tseq
        alpha.new <- 1/2 + sqrt(1 + 4*alpha*alpha)/2
        mom <- (alpha - 1)/alpha.new
        alpha <- alpha.new
        
        beta.vecy <- beta.new + mom*(beta.new - beta.old)
        objfn.val[k+1] <- objfn(beta.new, ...)
        beta.old <- beta.new
    }
    objfn.val <- objfn.val[!is.na(objfn.val)]
    value.obj <- objfn.val[k]
    return(list(par=beta.new, value.obj=value.obj, objfn.val=objfn.val, num.iter=k))
}
