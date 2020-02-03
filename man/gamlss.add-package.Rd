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
Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}. Cambridge. 

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.com/}).


Therneau T. M., Atkinson E. J. (2015) An Introduction to Recursive Partitioning Using the RPART Routines. Vignette in package rpart.

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer. 


Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press. 
}

\keyword{package}
\seealso{\code{\link[gamlss]{gamlss}}, \code{\link[gamlss.dist]{gamlss.family}}}  
\examples{
library(gamlss)
gn <- gamlss(R~ga(~te(Fl,A)), data=rent, family=GA) 
}
