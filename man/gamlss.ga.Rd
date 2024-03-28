\name{gamlss.ga}
\Rdversion{1.1}
\alias{gamlss.ga}
\alias{gamlss.ba}


\title{Support for Function ga() and ba()}
\description{
 This is support for the  smoother functions \code{ga()} and \code{ba()} intefaces for Simon Woood's \code{gam()} and \code{bam()}  functions from package \pkg{mgcv}.
It is not intended to be called directly by users. 
}
\usage{
gamlss.ga(x, y, w, xeval = NULL, ...)
gamlss.ba(x, y, w, xeval = NULL, ...)
}

\arguments{
  \item{x}{
the explanatory variables
}
  \item{y}{
iterative y variable
}
  \item{w}{
iterative weights
}
  \item{xeval}{
if xeval=TRUE then predicion is used
}
  \item{\dots}{
for extra arguments
}
}


\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, \bold{23}(7), 1--46, \doi{10.18637/jss.v023.i07}

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press. 
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby
}







% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
