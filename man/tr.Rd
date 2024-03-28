\name{tr}
\alias{tr}
\alias{gamlss.tr}

\title{
A interface function to use rpart() function within GAMLSS
}
\description{
The tr() function is a additive function to be used for GAMLSS models. It is an interface for the \code{rpart()} function of package \code{rpart}. The function tr() allows the user to use regression trees within gamlss. The great advantage of course comes from the fact GAMLSS models provide a variety of distributions and diagnostics. Note that the function gamlss.tr is not used by the user but it needed for the backfitting.
}
\usage{
tr(formula, method = c("rpart"), control = rpart.control(...), ...)
gamlss.tr(x, y, w, xeval = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A formula containing the expolanatory variables i.e. \code{~x1+x2+x3}.}
  \item{method}{only method "rpart" is supported at the moment}
  \item{control}{ control here is equivalent to \code{rpart.control()} function od package \code{rpart} }
  \item{x}{object passing informatio to the function}
  \item{y}{the iterative y variable}
  \item{w}{the iterative weights}
  \item{xeval}{whether prediction or not is used}
  \item{\dots}{additional arguments}
}

\details{
Note that, the gamlss fit maybe would not coverged. Also occasianly  the \code{gd.tol} argument in \code{gamlss} has to be increased. The 
}
\value{
Note that \code{tr} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function \code{additive.fit()} for backfitting which in turn uses \code{gamlss.tr()}
The result is a \code{rpart} object.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Ripley, B. D. (1996) \emph{Pattern Recognition and Neural Networks}. Cambridge. 


Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, \bold{23}(7), 1--46, \doi{10.18637/jss.v023.i07}

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{https://www.gamlss.com/}).

Therneau T. M., Atkinson E. J. (2015) An Introduction to Recursive Partitioning Using the RPART Routines. Vignette in package rpart.

Venables, W. N. and Ripley, B. D. (2002) \emph{Modern Applied Statistics with S}. Fourth edition. Springer. 
}
\author{
Mikis Stasinopoulos \email{mikis.stasinopoulos@gamlss.org}, Bob Rigby  based on work of Therneau and Atkison (2015)
}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{See Also as \code{\link{nn}}}
\examples{
data(rent)
#--- fitting gamlss+tree Nornal
library(rpart)
data(rent)
rg1 <- gamlss(R ~ tr(~A+Fl), data=rent, family=NO)
plot(rg1)
plot(getSmo(rg1))
text(getSmo(rg1))
\dontrun{
# fitting Gamma  errors
rg2 <- gamlss(R ~ tr(~A+Fl), data=rent, family=GA)
plot(rg2)
plot(getSmo(rg2))
text(getSmo(rg2))
#--- fitting also model in the variance
rg3 <- gamlss(R ~ tr(~A+Fl), sigma.fo=~tr(~Fl+A), data=rent, 
                family=GA, gd.tol=100, c.crit=0.1)
plot(rg3)
plot(getSmo(rg3))
text(getSmo(rg3))
plot(getSmo(rg3, what="sigma"))
text(getSmo(rg3, what="sigma"))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
