\name{blag}
\alias{blag}
\alias{llag}
\alias{wlag}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Functions to create lag values
}
\description{
The function \code{blag()} creates a basis for lag values of x, (a matrix of lag values of x). 
The function \code{llag()} creates a list with two components i) a basis matrix and  ii) weights to be used as prior weights in any regression analysis. The function \code{wlag()} can take a "mlags" object (created by \code{blag()}) or a vector and returns a vector with ones and zeros. This can be used as prior weights in any analysis which uses  \code{blag()}.  
}
\usage{
blag(x, lags = 1, from.lag=0, omit.na = FALSE, 
    value = NA, ...)
llag(x, ...)
wlag(x, lags = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
For \code{blag()} and \code{llag()} x is the vector for creating lags. For \code{wlag()} x is an \code{mlags} object created by \code{blag()}.   
}
  \item{lags}{
how many lags are required
}
 \item{from.lag}{where the lags are starting from. The default values is zero which indicates that the \code{x} is also included as a first column. If you want \code{x} not to included in the matrix use \code{from.lag=1}}
  \item{omit.na}{if true the first "lag" rows of the resulting matrix are omitted}
  \item{value}{value : what values should be set in the beginning of the lags columns, by default is set to NA}
  \item{\dots}{additional arguments}
}
\details{
Those three functions are design for helping a user to fit regression model using lags by generating the appropriate structures.
The function \code{blag()} creates a basis for lag values of x. It assumed that time runs from the oldest to the newest observations. That is, the latest observations are the most recent ones. The function \code{wlag()} take a  basis matrix of lags and creates a vector of weights which can be used as a prior weights for any regression type analysis which has the matrix as explanatory variable.  The function \code{llag()} creates a list with the matrix base for lags and the appropriate weights.   
}
\value{
The function \code{blag()} returns a "mlags" object (matrix of lag values).
The function \code{llag()} returns a list with components:
 \item{matrix}{The basis of the lag matrix}
 \item{weights}{The weights vector}
The function \code{wlag()} returns a vector of prior weights having, The vector starts with  zeros (as many as the number of lags)  and continues with ones.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos <\email{mikis.stasinopoulos@gamlss.org}>, Bob Rigby <\email{r.rigby@londonmet.ac.uk}> 
Vlasios Voudouris <\email{v.voudouris@londonmet.ac.uk}>, Majid  Djennad, Paul Eilers.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{  \code{\link{penLags}}}

\examples{
library(stats)
y <- arima.sim(500, model=list(ar=c(.4,.3,.1))) 
X <- blag(y, lags=5, from.lag=1,  value=0)
head(X)
w<-wlag(X)
library(gamlss)
m1<-gamlss(y~X, weights=w )
summary(m1)
plot(y)
lines(fitted(m1)~as.numeric(time(y)), col="blue")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
\keyword{ts}
