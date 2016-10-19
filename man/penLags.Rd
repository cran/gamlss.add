\name{penLags}
\alias{penLags}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Penalised Lag Regression Function}
\description{
The function \code{penLags()} fits a regression model to lags of an explanatory variable \code{x} or to lags of \code{y} itself.
The estimated coefficients of the lags are penalised using a quadratic penalty similar to P-splines.
}
\usage{
penLags(y, x, lags = 10, from.lag=0, weights = NULL, data = NULL, df = NULL, 
        lambda = NULL, start.lambda = 10, order = 1, 
        plot = FALSE, method = c("ML", "GAIC"), k = 2, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{The response variable}
  \item{x}{The explanatory variable which can be the response itself if autoregressive model is required}
  \item{lags}{The number of lags required}
 \item{from.lag}{from which lag value to start, the default is zero which means include the original \code{x} in the basis}
  \item{weights}{The prior weights}
  \item{data}{The data frame if needed}
  \item{df}{If not \code{NULL} this argument sets the required effective degrees of freedom for the penalty}
  \item{lambda}{If not \code{NULL} this argument sets the required smoothing parameter of  the penalty}
  \item{start.lambda}{Staring values for lambda for the local ML estimation}
  \item{order}{The order of the penalties in the beta coefficients}
  \item{plot}{Whether to plot the data and the fitted values}
  \item{method}{The method of estimating the smoothing parameter with two alternatives, i) \code{ML}: the local maximum likelihood estimation method (or PQL method) ii) \code{GAIC}: the generalised Akaike criterion method of estimating the smoothing parameter}
  \item{k}{The penalty required if the method  \code{GAIC} is used i.e.  \code{k=2} for AIC or \code{k=log(n)} if BIC (or SBC).}
  \item{\dots}{
for further arguments
}
}
\details{
This function is designed for fitting a simple  penalised lag regression model to a response variable. The meaning of simple in this case  is that only one explanatory variable can used (whether it is a true explanatory or the response variable itself) and only a normal assumption for the response is made.  For multiple explanatory variables and for different distributions within \code{gamlss} use the additive function \code{la}.  
}
\value{
Returns  \code{penLags} objects which has several method.   
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
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
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby,  Vlasios Voudouris,  Majid  Djennad, and Paul Eilers. 
}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

# generating data
y <- arima.sim(500, model=list(ar=c(.9,-.8)))
#----------------------------------
#fitting model with different order
m0 <- penLags(y,y, lag=20, order=0)
m1 <- penLags(y,y, lag=20, order=1)
m2 <- penLags(y,y, lag=20, order=2)
m3 <- penLags(y,y, lag=20, order=3)   
# chosing the order
AIC(m0, m1, m2, m3)
#---------------------------------
# look at the AR coefficients of the models
op <- par(mfrow=c(2,2))
 plot(coef(m0,"AR"), type="h")
 plot(coef(m1, "AR"), type="h")
 plot(coef(m2, "AR"), type="h")
 plot(coef(m3,"AR"), type="h")
par(op)
#-------------------------------
# refit and  plotting model
m1 <- penLags(y,y, lag=20, order=1, plot=TRUE)

# looking at the residuals
plot(resid(m1))
acf(resid(m1))
pacf(resid(m1))
# or better use plot, wp or dtop
plot(m1, ts=TRUE)
wp(m1)
dtop(m1)
# the coefficients
coef(m1)
coef(m1, "AR")
coef(m1, 'varComp')
#
print(m1)
#summary(m1)
# use prediction
plot(ts(c(y, predict(m1,100))))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
