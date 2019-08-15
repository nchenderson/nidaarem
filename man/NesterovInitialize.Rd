\name{NesterovInitialize}
\alias{NesterovInitialize}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Nesterov Initialization
}
\description{
Performs Nesterov initialization where an initial parameter value is found by
performing a preliminary number of Nesterov steps.
}
\usage{
NesterovInitialize(par, fixptfn, objfn = NULL, test = "monotone", ..., 
                   control = list(maxiter = 500))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{par}{
The initial parameter vector
}
  \item{fixptfn}{
The fixed point function which performs the basic parameter update.
}
  \item{objfn}{
The objective function of interest.
}
  \item{test}{
The test which determines termination of the Nesterov iterations. Can be
either "monotone" or "gradient".
}
  \item{\dots}{
Additional arguments that may be needed by \code{fixptfn} or \code{objfn}
}
  \item{control}{
A list of additional control parameters
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following components
  \item{par}{Initial value of the parameter vector}
  \item{value.obj}{value of objective function at termination}
  \item{objfn.val}{values of the objective function at each iteration}
  \item{num.iter}{number of fixed-point updates performed}
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Nicholas Henderson
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
LassoObj <- function(bbeta, X, y, lambda, stplngth) {
  X.beta <- as.vector(X\%*\%bbeta)
  ans <- sum((y - X.beta)*(y - X.beta))/2 + lambda*sum(abs(bbeta))
  return(-ans)
}

GDLassostep <- function(par, X, y, lambda, stplngth) {
  beta.vec <- par
  X.beta <- X\%*\%beta.vec
  beta.vec <- SoftThresh(beta.vec + stplngth*crossprod(X, y - X.beta), 
                         lambda=lambda*stplngth)
  return(beta.vec)
}

n <- 1000
p <- 100

beta.coef <- rt(p, df=2)
XX <- matrix(rnorm(n*p), nrow=n, ncol=p)
yy <- XX\%*\%beta.coef + rnorm(n) 
VV <- eigen(crossprod(XX))
LL <- max(VV$values)

ne.mon <- NesterovInitialize(par=rep(0, p), fixptfn=GDLassostep, objfn = LassoObj, 
                               test = "monotone", X=XX, y=yy, lambda=.0001, 
                               stplngth=1/LL)

ne.grad <- NesterovInitialize(par=rep(0, p), fixptfn=GDLassostep, objfn = LassoObj, 
                                test = "gradient", X=XX, y=yy, lambda=.0001, 
                                stplngth=1/LL)

ne.mon$num.iter
ne.grad$num.iter
}
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
