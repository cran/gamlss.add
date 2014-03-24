\name{tr}
\alias{tr}
\alias{gamlss.tr}

\title{
A interface function to use rpart() function within GAMLSS
}
\description{
The tr() function is a additive function to be used for GAMLSS models. It is an interface for the \code{rpart()} function of package \code{rpart}. The function tr() allows the user to use regression tres within gamlss. The great advantage of course comes from the fact GAMLSS models provide a variety of distributions and diagnostics. Note that the function gamlss.tr is not used by the user but it needed for the backfitting.
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
Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. 

}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk} based on work of Therneau and Atkison (19)
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
                family=GA, gd.tol=100, c,crit=0.1)
plot(rg3)
plot(getSmo(rg3))
text(getSmo(rg3))
plot(getSmo(rg3, what="sigma"))
text(getSmo(rg3, what="sigma"))}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}