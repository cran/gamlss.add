\name{la}
\alias{la}
\alias{la.control}
\alias{gamlss.la}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{
An Additive term for fitting Penalised Lags within a GAMLSS model
}
\description{
The function \code{la()} can be used for fitting penalised lags for  explanatory variables. 
}
\usage{
la(x, control = la.control(...), ...)
la.control(lags = 10, from.lag = 0, df = NULL, lambda = NULL, 
          start.lambda = 10, order = 1, plot = FALSE, 
         method = c("ML", "GAIC"), k = 2, ...)
gamlss.la(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the name of the explanatory variable}
  \item{y}{internally evaluated (iterative working variable)}
  \item{w}{internally evaluated (iterative weights)}
  \item{xeval}{internally evaluated no need for specification here}

  \item{control}{ a list of a number of control parameters for the fitting function \code{penLags()}}
  \item{lags}{the number of lags}
  \item{from.lag}{from which lag value to start, the default is zero which means include the original \code{x} in the basis}
  \item{df}{use this if you want to fix the degrees of freedom}
  \item{lambda}{use this if you would like to fix the smoothing parameter}
  \item{start.lambda}{initial starting value for lambda}
  \item{order}{the order of the penalty}
  \item{plot}{whether you would like a plot of the data}
  \item{method}{method of fitting if lambda or df are not specified}
  \item{k}{the penally used if method "GAIC" is used}
  \item{\dots}{for further arguments}
}
\details{
The idea of penalised lags is that we use a large amount of lags but we penalised their fitted coefficients and therefore we use few degrees of freedom.
The penally and method of fitting  are the same as in the \code{pb()} function of \code{gamlss}.  This function does not do the fitting this is achieved by the function \code{gamlss.la()} which uses the function \code{penLags} for the fitting}
\value{
a vector of zeros is returned, endowed with a number of attributes. The vector itself is used in the construction of the model matrix (contributing nothing), while the attributes are needed for the back-fitting algorithms of the additive fit.
}
\references{
Benjamin M. A., Rigby R. A. and Stasinopoulos D.M. (2003) Generalised Autoregressive Moving  Average Models.  \emph{J. Am. Statist. Ass.}, 98, 214-223.

Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}). 

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.}
\author{
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby,  Vlasios Voudouris, Majid  Djennad, and Paul Eilers. 
}

\note{
Note that an appropriate prior weight is needed in the gamlss fit
}


\seealso{\code{\link{penLags}}
}
\examples{
##  the data
 dax <- EuStockMarkets[,"DAX"]
 plot(dax)
## using a penalised autorgressive model
   w <- wlag(dax, lag=20)
  m1 <- gamlss(dax~ la(dax, lags=20, order=1, from.lag=1), weights=w) 
  lines(fitted(m1)~as.numeric(time(dax)), col=2)
 wp(m1, ylim.all=1) # not very good model
\dontrun{
## Try modelling the variance
  m2 <- gamlss(dax~ la(dax, lags=20, order=1, from.lag=1),
              sigma.fo=~la(dax^2, lags=10, order=1, from.lag=1), weights=w)
 wp(m2, ylim.all=1)# maybe the tails
  m3 <- gamlss(dax~ la(dax, lags=20, order=1, from.lag=1),
              sigma.fo=~la(dax^2, lags=10, order=1, from.lag=1), 
              weights=w,  family=TF)
  wp(m3, ylim.all=1) # better model 
  plot(m3, ts=TRUE)  # autocorrelation OK
## using FTSE to precict DAX           
ftse <- EuStockMarkets[,"FTSE"]
# fitting using penLags
  l1 <- penLags(dax, ftse, lags=30, plot=TRUE)
# similar model within gamlss  
   w <- wlag(ftse, lag=30)
  g1 <- gamlss(dax~ la(ftse, lags=30, order=1), weights=w)
lines(fitted(m1)~as.numeric(time(dax)))
  op <- par(mfrow=c(2,1))
# plotting the fitted coeficints of the AR terms  
plot(coef(l1, "AR"), type="h")
plot(coef(g1$mu.coefSmo[[1]])[-1], type="h")
par(op)
g2 <- gamlss(dax~ la(ftse, lags=30, order=1)+la(dax, lags=20, order=1, from.lag=1) , weights=w)
g3 <- gamlss(dax~ la(ftse, lags=30, order=1)+la(dax, lags=20, order=1, from.lag=1) ,
                sigma.fo=~la(dax^2, lags=10, order=1, from.lag=1), weights=w)
}}
\keyword{regression}
\keyword{ts}
