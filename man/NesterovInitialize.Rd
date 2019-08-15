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
NesterovInitialize(par, fixptfn, objfn = NULL, test = "monotone", ..., control = list(maxiter = 500))
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
  \item{par}{Description of 'comp1'}
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
a = 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
