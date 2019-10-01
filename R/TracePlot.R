TracePlot <- function(obj, opt.estimate=NULL, plot.trace=TRUE, col="black", lwd=2, ...) {
  if(is.null(opt.estimate)) {
    opt.estimate <- min(-obj$objfn.track)
  }
  tr <- -obj$objfn.track - opt.estimate
  niter <- length(obj$objfn.track)
  if(plot.trace) {
    plot(tr, type="n", log="y", las=1, xlab="Iteration", ylab="Objective - Optimum", ...)
    lines(1:niter, tr, col=col, lwd=lwd)
  }
  return(list(num.iter=niter, trace=tr))
}