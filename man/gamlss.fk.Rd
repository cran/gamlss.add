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

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{https://www.jstatsoft.org/article/view/v023i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).
}
\author{Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby}

\seealso{ \code{\link{fk}} }

\keyword{regression}
