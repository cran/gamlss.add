\name{gamlss.add-package}
\alias{gamlss.add-package}
\alias{gamlss.add}
\docType{package}
\title{
\packageTitle{gamlss.add}
}
\description{
\packageDescription{gamlss.add}
}
\details{

The DESCRIPTION file:
\packageDESCRIPTION{gamlss.add}
\packageIndices{gamlss.add}

}
\author{
\packageAuthor{gamlss.add}

Maintainer: \packageMaintainer{gamlss.add}
}
\references{
Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Therneau T. M., Atkinson E. J. (2015) An Introduction to Recursive Partitioning Using the RPART Routines. Vignette in package rpart.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. 


Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press. 
}

\keyword{package}
\seealso{\code{\link[gamlss]{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}}}  
\examples{
library(gamlss)
gn <- gamlss(R~ga(~te(Fl,A)), data=rent, family=GA) 
}
