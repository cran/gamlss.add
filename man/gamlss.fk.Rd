\name{gamlss.fk}
\Rdversion{1.1}
\alias{gamlss.fk}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Support for Function fk()
}
\description{
This is support for the functions \code{fk()}. 
It is not intended to be called directly by users. T
he function \code{gamlss.fk} is calling on the R function \code{curfit.free.knot()} of Sundar Dorai-Raj
}
\usage{
gamlss.fk(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the design matrix}
  \item{y}{the response variable}
  \item{w}{ prior weights}
  \item{xeval}{used in prediction } 
  \item{\dots}{for extra arguments}
}


\references{
Dierckx, P. (1991) \emph{Curve and Surface Fitting with Splines}, Oxford Science Publications

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} }

\seealso{ \code{\link{fk}} }

\keyword{regression}
