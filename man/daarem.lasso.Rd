\name{daarem.lasso}
\alias{daarem.lasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
NIDAAREM acceleration of penalized regression
}
\description{
Uses the NIDAAREM procedure to accelerate proximal gradient algorithms for fitting l1-penalized regression.
}
\usage{
daarem.lasso(par, X, y, lambda, stplngth = NULL, nesterov.init = FALSE, nesterov.safe = TRUE, family = c("gaussian", "binomial"), control = list())
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
initial value of the parameters
}
  \item{X}{
the n x p design matrix
}
  \item{y}{
vector of responses (should have length n)
}
  \item{lambda}{
the penalty term
}
  \item{stplngth}{
The step length used in the proximal gradient scheme.
}
  \item{nesterov.init}{
A logical indicator. If true, the NIDAAREM algorithm is used.
If false, the DAAREM algorithm is used.
}
  \item{nesterov.safe}{
Irrelevant - ignore for now
}
  \item{family}{
The link function used in the regression. Can be
either gaussian or binomial.
}
  \item{control}{
A list of additional control parameters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (par, X, y, lambda, stplngth = NULL, nesterov.init = FALSE, 
    nesterov.safe = TRUE, family = c("gaussian", "binomial"), 
    control = list()) 
{
    if ("objfn.check" \%in\% names(control)) {
        if (control$objfn.check) {
            control.default <- list(maxiter = 2000, order = 10, 
                tol = 1e-08, mon.tol = 1, cycl.mon.tol = 0, objfn.check = TRUE, 
                kappa = 25, alpha = 1.2)
        }
        else if (!control$objfn.check) {
            control.default <- list(maxiter = 2000, order = 10, 
                tol = 1e-08, mon.tol = 0.95, cycl.mon.tol = 0, 
                objfn.check = TRUE, kappa = 25, alpha = 1.2)
        }
    }
    else {
        control.default <- list(maxiter = 2000, order = 10, tol = 1e-08, 
            mon.tol = 1, cycl.mon.tol = 0, objfn.check = TRUE, 
            kappa = 25, alpha = 1.2)
    }
    namc <- names(control)
    if (!all(namc \%in\% names(control.default))) {
        stop("unknown names in control: ", namc[!(namc \%in\% names(control.default))])
    }
    control <- modifyList(control.default, control)
    family = match.arg(family)
    maxiter <- control$maxiter
    tol <- control$tol
    mon.tol <- control$mon.tol
    cycl.mon.tol <- control$cycl.mon.tol
    a1 <- control$alpha
    kappa <- control$kappa
    num.params <- length(par)
    nlag <- min(control$order, ceiling(num.params/2))
    lasso.pen <- lambda
    Fdiff <- Xdiff <- matrix(0, nrow = num.params, ncol = nlag)
    obj_funvals <- rep(NA, maxiter + 2)
    if (is.null(stplngth) & family == "gaussian") {
        stplngth <- FindStepLengthGaussian(X, y, lambda)
    }
    else if (is.null(stplngth) & family == "binomial") {
        stplngth <- FindStepLengthBinomial(X, y, lambda)
    }
    if (family == "gaussian" & control$objfn.check) {
        base_fn <- ifelse(n > p, "gauss_ngtp", "gauss_pgtn")
    }
    else if (family == "gaussian" & !control$objfn.check) {
        base_fn <- ifelse(n > p, "gauss_ngtp2", "gauss_pgtn2")
    }
    else if (family == "binomial" & control$objfn.check) {
        base_fn <- "binomial_b"
    }
    else if (family == "binomial" & !control$objfn.check) {
        base_fn <- "binomial_b2"
    }
    nest.fpevals <- 0
    if (nesterov.init & family == "gaussian" & nesterov.safe == 
        TRUE) {
        done <- FALSE
        count <- 0
        par.init <- par
        total.fpevals <- 0
        objfn.track <- NULL
        while (!done & count < 5) {
            neirun <- NesterovInitialize(par = par.init, fixptfn = GDLassoStep, 
                objfn = LassoObjFn, test = "monotone", X = X, 
                y = y, lambda = lambda, stplngth = stplngth, 
                control = list(maxiter = maxiter))
            par <- neirun$par
            nest.fpevals <- neirun$num.iter
            if (count == 0) {
                maxiter.nest <- max(5 * nest.fpevals, 1000)
            }
            ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, 
                X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                  nest.fpevals, tol, mon.tol, cycl.mon.tol), 
                gauss_pgtn = daarem_lasso_gaussian_pgtn(par, 
                  X, y, lambda, stplngth, nlag, a1, kappa, maxiter.nest, 
                  tol, mon.tol, cycl.mon.tol), gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, 
                  X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                    nest.fpevals, tol, mon.tol, cycl.mon.tol), 
                gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, 
                  X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                    nest.fpevals, tol, mon.tol, cycl.mon.tol), 
                binomial_b = daarem_lasso_binomial(par, X, y, 
                  lambda, stplngth, nlag, a1, kappa, maxiter - 
                    nest.fpevals, tol, mon.tol, cycl.mon.tol), 
                binomial_b2 = daarem_lasso_binomial2(par, X, 
                  y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                    nest.fpevals, tol, mon.tol, cycl.mon.tol))
            total.fpevals <- total.fpevals + ans$fpevals + nest.fpevals
            objfn.track <- c(objfn.track, neirun$objfn.track, 
                ans$objfn.track)
            done <- ans$convergence
            count <- count + 1
            par.init <- ans$par
        }
        ans$fpevals <- total.fpevals
        ans$objfn.track <- objfn.track
    }
    else if (nesterov.init & family == "gaussian" & nesterov.safe == 
        FALSE) {
        par.init <- par
        neirun <- NesterovInitialize(par = par.init, fixptfn = GDLassoStep, 
            objfn = LassoObjFn, test = "monotone", X = X, y = y, 
            lambda = lambda, stplngth = stplngth, control = list(maxiter = maxiter))
        nest.fpevals <- neirun$num.iter
        par <- neirun$par
        ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_pgtn = daarem_lasso_gaussian_pgtn(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), binomial_b = daarem_lasso_binomial(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), binomial_b2 = daarem_lasso_binomial2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol))
        ans$fpevals <- ans$fpevals + nest.fpevals
        ans$objfn.track <- c(neirun$objfn.track, ans$objfn.track)
    }
    else {
        ans <- switch(base_fn, gauss_ngtp = daarem_lasso_gaussian_ngtp(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_pgtn = daarem_lasso_gaussian_pgtn(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_ngtp2 = daarem_lasso_gaussian_ngtp2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), gauss_pgtn2 = daarem_lasso_gaussian_pgtn2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), binomial_b = daarem_lasso_binomial(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol), binomial_b2 = daarem_lasso_binomial2(par, 
            X, y, lambda, stplngth, nlag, a1, kappa, maxiter - 
                nest.fpevals, tol, mon.tol, cycl.mon.tol))
    }
    return(ans)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ acceleration }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
