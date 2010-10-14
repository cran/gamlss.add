\name{fk}
\Rdversion{1.1}
\alias{fk}

\title{A  function to fit break points within GAMLSS}

\description{
The \code{fk()} function is a additive function to be used for GAMLSS models. 
It is an interface for the \code{fitFreeKnots()} function of package 
\pkg{gamlss.util}. The functions \code{fitFreeKnots()} was first based on the  \code{curfit.free.knot()} function of package 
DierckxSpline of  Sundar Dorai-Raj and Spencer Graves. The function \code{fk()} allows the user to 
use the free knots function \code{fitFreeKnots()} within \code{gamlss}. 
The great advantage of course comes from the fact GAMLSS models provide a variety of distributions and diagnostics.
}
\usage{
fk(x, degree = 3, start = NULL, ...)
}


\arguments{
  \item{x}{the x-variable}
  \item{degree}{the degree of the spline function fitted}
  \item{start}{starting values for the breakpoints. If are set the number of break points is also determined
                by the length of \code{start}}
 \item{\dots}{for extra arguments}
}
\details{
Note that \code{fk} itself does no smoothing; it simply sets things up for the function \code{gamlss()} 
which in turn uses the function  \code{additive.fit()} for backfitting which in turn uses \code{gamlss.fk()}.
Note that, finding the break points is not a trivial problem and therefor  multiple maximum points can occur.  
More details about the free knot splines can be found in  Dierckx, (1991).

The \code{gamlss} algorithm used a modified backfitting in this case, that is, it fits the linear part fist.  
Note that trying to predict outside the x-range can be dangerous as the example below shows.  
}
\value{
The \code{gamlss} object saved contains the last fitted object which can be accessed using 
\code{obj$par.coefSmo} where \code{obj} is the fitted \code{gamlss} object \code{par} is the relevant distribution
 parameter.}
   
\references{

Dierckx, P. (1991) \emph{Curve and Surface Fitting with Splines}, Oxford Science Publications

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}
}


\seealso{
  \code{\link{gamlss.fk}}} 


\examples{
library(gamlss.util)
# creating  a linear + linear function
    x <- seq(0,10, length.out=201)
 knot <- 5
 set.seed(12543)
    mu <- ifelse(x<=knot,5+0.5*x,5+0.5*x+(x-knot))
    y <- rNO(201, mu=mu, sigma=.5)
# plot the data
 plot(y~x, xlim=c(-1,13), ylim=c(3,17))
# fit model using curfit
 m1 <- fitFreeKnots(x, y, knots=4, degree=1)
# fitted values
lines(fitted(m1)~x, col="red", lwd="3")
# predict
pm1<-predict(m1, newdata=-1:12)
 points(-1:12,pm1, col="red",pch = 21, bg="blue")
#------------------------------------------------
# now gamlss
#------------------------------------------------
# get the fk function
# fit model
g0 <- gamlss(y~fk(x, degree=1, start=c(4)))
# creating data frame
da <- data.frame(y,x)
g1 <- gamlss(y~fk(x, degree=1, start=c(4)), data=da)
nd<-data.frame(x=-1:12)
pg1<-predict(g1,  newdata=nd)
lines(fitted(g1)~x, col="purple", lwd="3")
points(pg1~I(-1:12), col="darkgreen", pch = 21, bg="green" )
# note the predicting values outside the x-range
# it predicts the global linear fit
#-------------------------------------------------
#-------------------------------------------------
# now negative binomial data 
eta1 <- ifelse(x<=knot,5+0.5*x,5+0.5*x+(x-knot))
 set.seed(143)
 y <- rNBI(201, mu=exp(eta1/7), sigma=.1)
 da <- data.frame(y=y,x=x)
#rm(y,x)
plot(y~x, data=da)
# fitting the model in gamlss using profile deviance
n1 <- quote(gamlss(y ~ x+I((x>this)*(x-this)),family=NBI,data=da))
prof.term(n1, min=1, max=9, step=.1, criterion="GD")
# now fit the model using fk
g1 <- gamlss(y~fk(x, degree=1, start=c(4)), data=da, family=NBI)
# get the breakpoint
knots(g1$mu.coefSmo[[1]])
lines(fitted(g1)~x, data=da)
#------------------------------------------------
# the aids data
# using fk()
data(aids)
a1<-gamlss(y~fk(x, degree=1, start=25)+qrt, data=aids, family=NBI)
knots(a1$mu.coefSmo[[1]])
# using profile deviance
aids.1 <- quote(gamlss(y ~ x+I((x>this)*(x-this))+qrt,family=NBI,data=aids))
prof.term(aids.1, min=16, max=21, step=.1, criterion="GD")
# using non-linear estimation
library(gamlss.nl)
naids <- nl.obj(~p1*(x-p2)*I((x>p2)), start=c(0.2, 20), data=aids)
mod1 <- gamlss(y~nl(naids)+x+qrt, data=aids , family=NBI) 
mod1$mu.coefSmo 
GAIC(mod1,a1) # almost identical
# plotting the fit
with(aids, plot(x,y,pch=21,bg=c("red","green3","blue","yellow")[unclass(qrt)]))
lines(fitted(a1)~aids$x)
#-------------------------------------------------
}

\keyword{regression}
