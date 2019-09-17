DampingFind <- function(uy.sq, dvec, aa, kappa, sk, Ftf, lambda.start=NULL,
                        r.start=NULL, maxit=10) {
  ## uy.sq is a vector whose jth component (U^t*y)^2
  ## d.sq is a vector whose jth component is d_j^2.
  ## (aa, kappa, sk) - parameters in determining delta_k
  ## Ftf <- F_{k}^T f_{k}
  ##
  ## Output: value of penalty term such that (approximately)
  ## .       || beta.ridge ||^2/||beta.ols||^2 = delta_{k}^{2}, with 0 < target < 1.

  if(is.null(lambda.start) | is.na(lambda.start)) {
    lambda.start <- 100
  }
  if(is.null(r.start) | is.na(r.start)) {
    r.start <- 0
  }
  if(sum(dvec==0.0) > 0) {
      ind <- dvec > 0
      dvec <- dvec[ind]
      uy.sq <- uy.sq[ind]
      ## what about if dvec has length 1?
  }
  pow <- kappa - sk
  target <- exp(-0.5*log1p(aa^pow))  ### This is sqrt(delta_k)
  #d.sq <- pmax(dvec*dvec, 1e-7)  ## in case the least-squares problem is very poorly conditioned
  d.sq <- dvec*dvec
  betahat.ls <- uy.sq/d.sq
  betahat.ls.norm <- sqrt(sum(betahat.ls))
  vk <- target*betahat.ls.norm

  if(vk == 0) {
      ## the norm of the betas is zero, so the value of lambda shouldn't matter
      return(list(lambda=lambda.start, rr=r.start))
  }
  ### Initialize lambda and lower and upper bounds
  lambda <- lambda.start - r.start/vk
  LL <- (betahat.ls.norm*(betahat.ls.norm - vk))/sum(uy.sq/(d.sq*d.sq))
  UU <- Ftf/vk

  ## Compute lstop and ustop
  pow.low <- pow + 0.5
  pow.up <- pow - 0.5
  l.stop <- exp(-0.5*log1p(aa^pow.low))
  u.stop <- exp(-0.5*log1p(aa^pow.up))

  ### Start iterations
  for(k in 1:maxit) {
    ### Check if lambda is inside lower and upper bounds
    if(lambda <= LL | lambda >= UU) {
      lambda <- max(.0001*UU, sqrt(LL*UU))
    }
    ### Evaluate ||s(lambda)|| and \phi(lambda)/phi'(lambda)

    d.lambda <- (dvec/(d.sq + lambda))^2
    d.prime <- d.lambda/(d.sq + lambda)
    s.norm <- sqrt(sum(uy.sq*d.lambda))
    phi.val <- s.norm - vk
    phi.der <- (-1)*sum(uy.sq*d.prime)/s.norm
    phi.ratio <- phi.val/phi.der
   
    d.u <- (dvec/(d.sq + LL))^2
    s.up <- sqrt(sum(uy.sq*d.u))
    ### Initial Lower bound is not correct

    ### Check convergence
    if(s.norm <= u.stop*betahat.ls.norm & s.norm >= l.stop*betahat.ls.norm) {
      break
    }
    
    ### If not converged, update lower and upper bounds
    UU <- ifelse(phi.val >= 0, UU, lambda)
    LL <- max(LL, lambda - phi.ratio)

    ### Now update lambda
    lambda <- lambda - (s.norm*phi.ratio)/vk
    bb <- (s.norm*phi.ratio)/vk
  }
  ans <- list(lambda=lambda, rr=s.norm*phi.ratio)
  return(ans)
}
