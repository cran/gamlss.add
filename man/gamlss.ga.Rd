\name{gamlss.ga}
\Rdversion{1.1}
\alias{gamlss.ga}

\title{
Support for Function ga() 
}
\description{
 This is support for the  smoother function ga() an inteface for Simon Woood's \code{gam()} function.
It is not intended to be called directly by users. 
}
\usage{
gamlss.ga(x, y, w, xeval = NULL, ...)
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

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press. 
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}
}







% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
