\name{centilesTwo}
\alias{centilesTwo}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Centiles contour plots in GAMLSS}
\description{
This function \code{centilesTwo()} plots two dimensional centiles contour plots for GAMLSS models. 
}
\usage{
centilesTwo(object, grid.x1, grid.x2, x1.name, x2.name, 
       cent = 0.05, dist = 0.01, points = TRUE, 
       other = list(), point.col = 1, point.pch = ".", 
       image = FALSE, image.col = heat.colors(12), ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{an gamlss object}
  \item{grid.x1}{grid values for x-variable one}
  \item{grid.x2}{grid values for x-variable two}
  \item{x1.name}{the name of x-variable on}
  \item{x2.name}{the name of x-variable two}
  \item{cent}{the required centiles}
  \item{dist}{the distance}
  \item{points}{whether to plot the data points}
  \item{other}{a list having other explanatory variables at fixed values}
  \item{point.col}{the colour of the data points}
  \item{point.pch}{the type of the data point}
  \item{image}{whether to plot using the \code{image9} function}
  \item{image.col}{the colour scheme }
  \item{\dots}{for extra arguments for the \code{contour()} function}
}
\details{
The function uses the function \code{exclude.too.far()} of the package \pkg{mgcv}.
}
\value{
Produce a contour plot.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Rigby R.A., Stasinopoulos D. M., Heller G., and De Bastiani F., (2019) \emph{Distributions for Modeling Location, Scale and Shape: Using GAMLSS in R}, Chapman and Hall/CRC.

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017) \emph{Flexible Regression and Smoothing: Using GAMLSS in R}, Chapman and Hall/CRC. 

(see also \url{http://www.gamlss.com/}).

Wood S.N. (2006) \emph{Generalized Additive Models: An Introduction with R}. Chapman and Hall/CRC Press. 
}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby, Fernanda De Bastiani
}

\seealso{
\code{\link{centiles}}
}
\examples{
\dontrun{ 
data(plasma)
m1 <- gamlss(betadiet ~ ga(~te(age, fiber)), sigma.formula = ~1,  
    nu.formula = ~ga(~te(age, fiber)), tau.formula = ~1,  
    family = BCTo, data = plasma)
centilesTwo(m1, 18:90, seq(2.5,38, 0.5), age, fiber, cent=0.05, dist=.1,
            xlab="age", ylab='fiber') 
centilesTwo(m1, 18:90, seq(2.5,38, 0.5), age, fiber, cent=0.95, dist=.1)
}}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{regression}
