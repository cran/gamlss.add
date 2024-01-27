\name{fitFixedKnots}

\alias{fitFixedKnots}
\alias{fitFreeKnots}
\title{Functions to Fit Univariate Break Point Regression Models}
\description{
There are two main functions here. The functions \code{fitFixedKnots} allows the fit a univariate regression using piecewise polynomials with "known" break points while the function \code{fitFreeKnots} estimates the break points.  
}
\usage{
fitFixedKnots(y, x, weights = NULL, knots = NULL, data = NULL, degree = 3, 
             fixed = NULL, base=c("trun","Bbase"), ...)
fitFreeKnots(y, x, weights = NULL, knots = NULL, degree = 3, fixed =
                 NULL, trace = 0, data = NULL, base=c("trun","Bbase"), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{the x variable (explanatory)}
  \item{y}{the response variable}
  \item{weights}{the prior weights}
  \item{knots}{the position of the interior knots for  \code{fitFixedKnots} or starting values for \code{fitFreeKnots}}
  \item{data}{the data frame}
  \item{degree}{the degree if the piecewise  polynomials}
  \item{fixed}{this is to be able to fit fixed break points}
  \item{base}{The basis for the piecewise polynomials, \code{turn} for truncated (default) and \code{Bbase} for B-base piecewise polynomials   }
  \item{trace}{controlling the trace of of \code{optim()} }
  \item{\dots}{for extra arguments}
  }
\details{
The functions \code{fitFreeKnots()} is loosely based on the  \code{curfit.free.knot()} function of package 
\pkg{DierckxSpline} of  Sundar Dorai-Raj and Spencer Graves. 
}
\value{
The functions \code{fitFixedKnots} and \code{fitFreeKnots} return an object \code{FixBreakPointsReg} and 
\code{FreeBreakPointsReg} respectively with the following items:
 \item{fitted.values}{the fitted values of the model}
 \item{residuals}{the residuals of the model}
 \item{df}{the degrees of freedom fitted in the model}
 \item{rss}{the residuals sum of squares}
 \item{knots}{the knots used in creating the beta-function base}
 \item{fixed}{the fixed break points if any}
 \item{breakPoints}{the interior (estimated) break points (or knots)}
 \item{coef}{the coefficients of the linear  part of the model}
 \item{degree}{the degree of the piecewise polynomial}
 \item{y}{the y variable}
 \item{x}{the x variable}
 \item{w}{the prior weights}
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
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}
}
\note{
The prediction function in piecewise polynomials using the B-spline basis is tricky because by adding the newdata for x to the current one the B-basis function for the piecewise polynomials changes.  This does not seems to be the case with the truncated basis, that is, when the option \code{base="turn"} is used (see the example below).

If the newdata are outside the range of the old x then there could a considerable discrepancies between the all fitted values and the predicted ones if the option \code{base="Bbase"} is used. The prediction function for the objects \code{FixBreakPointsReg} or \code{FreeBreakPointsReg}
has the option \code{old.x.range=TRUE} which allow the user two choices:  

The first is to use the old end-points for the creation of the new B-basis which were determine from the original range of x. This choice is implemented as a default in the \code{predict} method for \code{FixBreakPointsReg} and \code{FreeBreakPointsReg} objects with the argument \code{old.x.range=TRUE}. 

The second  is to create new end-points from the new and old data  x values. In this case 
the range of x  will be bigger that the original one if the \code{newdata} has values outside the original x range.  
In this case (\code{old.x.range=FALSE}) the prediction could be possible  better outside the x range but would not coincide with the original predictions i.e. fitted(model)
since basis have changed. 



}

%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{
# creating  a linear + linear function
   x <- seq(0,10, length.out=201)
knot <- 5
 set.seed(12543)
 mu <- ifelse(x<=knot,5+0.5*x,5+0.5*x+(x-knot))
  y <- rNO(201, mu=mu, sigma=.5)
# plot the data
 plot(y~x, xlim=c(-1,13), ylim=c(3,18))

# fit model using fixed break points
 m1 <- fitFixedKnots(y, x, knots=5, degree=1)
knots(m1)
lines(fitted(m1)~x, col="red")

# now estimating the knot
m2 <- fitFreeKnots(y, x, knots=5, degree=1)
knots(m2)
summary(m2)

# now predicting 
plot(y~x, xlim=c(-5,13), ylim=c(3,18))
lines(fitted(m2)~x, col="green", lwd=3)
points(-2:13,predict(m2, newdata=-2:13), col="red",pch = 21, bg="blue")
points(-2:13,predict(m2, newdata=-2:13, old.x.range=FALSE), col="red",pch = 21, bg="grey")

# fit different basis 
m21 <- fitFreeKnots(y, x, knots=5, degree=1, base="Bbase")
deviance(m2) 
deviance(m21) # should be identical

# predicting with m21 
 plot(y~x, xlim=c(-5,13), ylim=c(3,18))
lines(fitted(m21)~x, col="green", lwd=3)
points(-2:13,predict(m21, newdata=-2:13), col="red",pch = 21, bg="blue")
points(-2:13,predict(m21, newdata=-2:13, old.x.range=FALSE), col="red",pch = 21, bg="grey")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
