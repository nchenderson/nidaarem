\name{nidaarem}
\alias{nidaarem}
\title{
Nesterov Initialized DAAREM
}
\description{
Performs the DAAREM acceleration scheme with Nesterov
initialization and alternating monotonicity
}
\usage{
nidaarem(par, fixptfn, objfn, ..., control = list())
}
\arguments{
  \item{par}{
A vector of starting values.
}
  \item{fixptfn}{
A vector function that performs the fixed-point mapping.
}
  \item{objfn}{
A scalar function that evaluates the objective function to maximized.
}
  \item{\dots}{
Arguments passed to \code{fixptfn} and \code{objfn}
}
  \item{control}{
A list of control parameters specifying any changes to default values of algorithm
control parameters. Full names of control list elements must be specified,
otherwise, user-specifications are ignored. See *Details*.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
A list with the following components:
\item{par}{Parameter, \eqn{x*}{x*} that are the fixed-point of \eqn{G}{G}
  such that \eqn{x* = G(x*)}{x* = G(x*)}, if convergence is successful.}

\item{value.objfn}{The value of the objective function $L$ at termination.}

\item{fpevals}{Number of times the fixed-point function \code{fixptfn} was evaluated.}

\item{objfevals}{Number of times the objective function \code{objfn} was evaluated.}

\item{convergence}{An integer code indicating type of convergence.  \code{0}
  indicates successful convergence, whereas \code{1} denotes failure to
  converge.}
\item{objfn.track}{A vector containing the value of the objective function at each iteration.}
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

PP <- matrix(rnorm(9), nrow=3, ncol=3)
PP[c(3,4,7,8)] <- NA
ind <- which(!is.na(PP))

A <- PP[ind]
PA <- PP
PA[-ind] <- 0.0

par0 <- rnorm(9)

matcompl.test <- nidaarem(par=par0, fixptfn=MatrixCompleteUpdate, objfn=MatrixCompleteObj, 
                          lambda=2, A=A, PA=PA, ind=ind)

}
\keyword{ acceleration }% use one of  RShowDoc("KEYWORDS")
\keyword{ optimization }% __ONLY ONE__ keyword per line
