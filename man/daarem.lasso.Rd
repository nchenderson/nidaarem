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
daarem.lasso(par, X, y, lambda, stplngth = NULL, nesterov.init = FALSE, family = c("gaussian", "binomial"), control = list())
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
a = 2
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ acceleration }% use one of  RShowDoc("KEYWORDS")
\keyword{ regression }% __ONLY ONE__ keyword per line
