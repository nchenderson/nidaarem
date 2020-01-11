\name{Nesterov}
\alias{Nesterov}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nesterov Acceleration with Restarts
}
\description{
Performs acceleration of an arbitrary fixed-point update using Nesterov's momentum method.
An option to use Nesterov momentum with restarts is included.
}
\usage{
Nesterov(par, fixptfn, objfn = NULL, restart = FALSE, ..., control = list(maxiter = 1000, tol = 1e-07))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
Initial parameter vector.
}
  \item{fixptfn}{
Function performing fixed-point update.
}
  \item{objfn}{
Objective function that one wants to maximize.
}
  \item{restart}{
A logical. Indicates whether or not restarts should be used in the acceleration scheme.
}
  \item{\dots}{
Additional arguments that could be passed to the fixed-point function.
}
  \item{control}{
A list of control parameters.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
alist(par=beta.old, value.objfn=value.obj, objfn.track=objfn.val, fpevals=k, convergence=conv)
A list with the following components:
   \item{par}{the final parameter values}
   \item{value.objfn}{the final value of the objective function}
   \item{objfn.track}{the value of the objective function at each iteration}
   \item{fpevals}{the number of evaluations of the fixed-point function}
   \item{convergence}{An indicator of whether or not convergence was achieved.}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson and Ravi Varadhan
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
n <- 1000
p <- 100

beta.coef <- rt(p, df=2)
XX <- matrix(rnorm(n*p), nrow=n, ncol=p)
yy <- XX\%*\%beta.coef + rnorm(n) 
VV <- eigen(crossprod(XX))
LL <- max(VV$values)

ne.test <- Nesterov(par=rep(0, p), fixptfn=GDLassoStep, objfn = LassoObjFn, 
                    restart=TRUE, X=XX, y=yy, lambda=.0001, stplngth=1/LL)
ne.test$fpevals
TracePlot(ne.test)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ acceleration }% use one of  RShowDoc("KEYWORDS")
\keyword{ restarts }% __ONLY ONE__ keyword per line
