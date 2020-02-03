\name{gamlss.nn}
\Rdversion{1.1}
\alias{gamlss.nn}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Support for Function nn() 
}
\description{
This is support for the  smoother function nn() an interface for Brian Reply's \code{nnet()} function.
It is not intended to be called directly by users. 
}
\usage{
gamlss.nn(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the explanatory variables}
  \item{y}{iterative y variable}
  \item{w}{iterative weights}
  \item{xeval}{if xeval=TRUE then predicion is used}
  \item{\dots}{for extra arguments}
}

\references{

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}. Cambridge. 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.com/}).

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer. 


}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby
}


\seealso{
\code{\link{fk}}}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
