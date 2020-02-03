\name{ga}
\Rdversion{1.1}
\alias{ga}
\alias{ga.control}
\alias{ba}
\alias{ba.control}

\title{A interface functions to use Simon Wood's gam() and bam() functions within GAMLSS}
\description{
The \code{ga()} and \code{ba()}  functions are a additive functions to be used within GAMLSS models. 
They are  interfaces for the \code{gam()} and the \code{bam()} functions of package  \code{mgcv} of Simon Wood. The functions  \code{gam()} and the \code{bam()} allows the user to use all the available smoothers of the package \code{mcgv()} within \code{gamlss}. The great advantage of course come from fitting models outside the 
exponential family.    
}
\usage{
ga(formula, control = ga.control(...), ...)

ba(formula, control = ba.control(...), ...)

ga.control(offset = NULL, method = "REML",
         optimizer = c("outer", "newton"), control = list(), 
         scale = 0, select = FALSE, knots = NULL, 
         sp = NULL, min.sp = NULL, H = NULL, gamma = 1, 
         paraPen = NULL, in.out = NULL, 
         drop.unused.levels = TRUE, drop.intercept = NULL, 
         discrete = FALSE, ...)
     
ba.control(offset = NULL, method = "fREML", control = list(), 
          select = FALSE, scale = 0, gamma = 1, knots = NULL, 
          sp = NULL, min.sp = NULL, paraPen = NULL, 
          chunk.size = 10000, rho = 0, AR.start = NULL, 
          discrete = TRUE, cluster = NULL, nthreads = 2, 
          gc.level = 1, use.chol = FALSE, samfrac = 1, 
          coef = NULL, drop.unused.levels = TRUE, 
          drop.intercept = NULL, ...)     

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ A formula containing \code{s()} and \code{te} functions i.e. ~s(x1)+ te(x2,x3).}
  \item{offset}{the offset in the formula}
  \item{method}{the \code{method} argument in \code{gam()} and \code{bam()}}
  \item{optimizer}{the method optimizer in \code{gam()}}
  \item{control}{values for the \code{gam.control()}}
  \item{scale}{for the scale parameter}
  \item{select}{the \code{select} argument in \code{gam()} and \code{bam()}}
  \item{knots}{the \code{knots} argument in \code{gam()} and \code{bam()}}  
  \item{sp}{the \code{sp} argument in \code{gam()} and \code{bam()}}
  \item{min.sp}{the \code{min.sp} argument in \code{gam()} and \code{bam()}}
  \item{H}{a user supplied fixed quadratic penalty on the parameters in \code{gam()}}
  \item{gamma}{the \code{gamma} argument in \code{gam()} and \code{bam()}}
  \item{paraPen}{the \code{paraPen} argument in \code{gam()} and \code{bam()}}
  \item{in.out}{the \code{in.out} argument in \code{gam()}}
  \item{drop.unused.levels}{by default unused levels are dropped from factors before fitting for \code{gam()} and \code{bam()}}
  \item{drop.intercept}{set to TRUE to force the model to really not have the a constant in the parametric model part for \code{gam()} and \code{bam()}}
  \item{discrete}{see \code{bam} and \code{gam} for details}
  \item{chunk.size}{see the help for \code{bam()}.}
  \item{rho}{for an \code{AR1} error model, see the help for \code{bam()}  }
  \item{AR.start}{for an \code{AR1} error model, see the help for \code{bam()}}
  \item{cluster}{see the help for \code{bam()}}
  \item{nthreads}{
Number of threads to use for non-cluster computation see  the help for \code{bam()}}
  \item{gc.level}{keepingf the memory footprint down, see  the help for \code{bam()}}
  \item{use.chol}{see  the help for \code{bam()}}
  \item{samfrac}{
see  the help for \code{bam()}
}
  \item{coef}{initial values for model coefficients}
  \item{\dots}{extra options to pass to gam.control()}
}
\details{
Note that \code{ga} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function 
  \code{additive.fit()} for back-fitting which in turn uses \code{gamlss.ga()}

Note that, in our (limited) experience,  for normal errors or exponential family, the  fitted models using \code{gam()} 
and \code{ga()} within \code{gamlss()} are identical
or at least very similar.  This is particularly true if the default values for \code{gam()} are used.   
   
}
\value{ the fitted values of the smoother is returned, endowed with a number of attributes. 
The smoother fitted values are used in the construction of the overall fitted values of the particular distribution parameter.
The attributes can be use to obtain information about the individual fit. In particular the \code{coefSmo} within the parameters
of the fitted model contains the final additive fit.
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
\author{Mikis Stasinopoulos, \email{d.stasinopoulos@londonmet.ac.uk} }

\section{Warning }{The function is experimental so please report any peculiar behaviour to the authors} 


\examples{
library(mgcv)
data(rent)
#---------------------------------------------------------
## normal errors one x-variable
ga1 <- gam(R~s(Fl, bs="ps", k=20), data=rent, method="REML")
gn1 <- gamlss(R~ga(~s(Fl, bs="ps", k=20), method="REML"), data=rent) # additive
gb1 <- gamlss(R~pb(Fl), data=rent) # additive
AIC(ga1,gn1, gb1, k=0)
AIC(ga1,gn1, gb1)
#--------------------------------------------------------
## normal error additive in Fl and A
ga2 <- gam(R~s(Fl)+s(A), method="REML", data=rent)
gn2 <- gamlss(R~ga(~s(Fl)+s(A), method="REML"),  data=rent) # additive
gb2 <- gamlss(R~pb(Fl)+pb(A), data=rent) # additive
AIC(ga2,gn2, gb2, k=0)
AIC(ga2,gn2, gb2)
#---------------------------------------------------------
\dontrun{
## gamma error additive in Fl and A
ga3 <- gam(R~s(Fl)+s(A), method="REML", data=rent, family=Gamma(log))
gn3 <- gamlss(R~ga(~s(Fl)+s(A), method="REML"), data=rent, family=GA)# additive
gb3 <- gamlss(R~pb(Fl)+pb(A), data=rent, family=GA) # additive
AIC(ga3,gn3, gb3, k=0)
AIC(ga3,gn3, gb3)
#---------------------------------------------------------
## gamma error surface fitting
ga4 <-gam(R~s(Fl,A), method="REML", data=rent, family=Gamma(log))
gn4 <- gamlss(R~ga(~s(Fl,A), method="REML"), data=rent, family=GA) 
AIC(ga4,gn4, k=0)
AIC(ga4,gn4)
## plot the fitted surfaces
op<-par(mfrow=c(1,2))
vis.gam(ga4)
vis.gam(getSmo(gn4))
par(op)
## contour plot using mgcv's plot() function
plot(getSmo(gn4))
#---------------------------------------------------------
## predict
newrent <- data.frame(expand.grid(Fl=seq(30,120,5), A=seq(1890,1990,5 )))
newrent1 <-newrent2 <- newrent
newrent1$pred <- predict(ga4, newdata=newrent, type="response")
newrent2$pred <- predict(gn4, newdata=newrent, type="response")
library(lattice)
wf1<-wireframe(pred~Fl*A, newrent1, aspect=c(1,0.5), drape=TRUE, 
             colorkey=(list(space="right", height=0.6)), main="gam()")
wf2<-wireframe(pred~Fl*A, newrent2, aspect=c(1,0.5), drape=TRUE, 
          colorkey=(list(space="right", height=0.6)), main="gamlss()")
print(wf1, split=c(1,1,2,1), more=TRUE)
print(wf2, split=c(2,1,2,1))
#---------------------------------------------------------
##gamma error two variables te() function
ga5 <-  gam(R~te(Fl,A), data=rent, family=Gamma(log))
gn5 <- gamlss(R~ga(~te(Fl,A)), data=rent, family=GA) 
AIC(ga5,gn5)
AIC(ga5,gn5, k=0)
op<-par(mfrow=c(1,2))
vis.gam(ga5)
vis.gam(getSmo(gn5))
par(op)
#----------------------------------------------------------
## use of Markov random fields 
## example from package mgcv of Simon Wood
## Load Columbus Ohio crime data (see ?columbus for details and credits)
data(columb)       ## data frame
data(columb.polys) ## district shapes list
xt <- list(polys=columb.polys) ## neighbourhood structure info for MRF
## First a full rank MRF...
b <- gam(crime ~ s(district,bs="mrf",xt=xt),data=columb,method="REML")
bb <- gamlss(crime~ ga(~s(district,bs="mrf",xt=xt), method="REML"), data=columb)
AIC(b,bb, k=0)
op<-par(mfrow=c(2,2))
plot(b,scheme=1)
plot(bb$mu.coefSmo[[1]], scheme=1)
## Compare to reduced rank version...
b <- gam(crime ~ s(district,bs="mrf",k=20,xt=xt),data=columb,method="REML")
bb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt), method="REML"), 
             data=columb)
AIC(b,bb, k=0)
plot(b,scheme=1)
plot(bb$mu.coefSmo[[1]], scheme=1)
par(op)
## An important covariate added...
b <- gam(crime ~ s(district,bs="mrf",k=20,xt=xt)+s(income),
         data=columb,method="REML")
## x in gam() 
bb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt)+s(income), 
             method="REML"), data=columb)
## x in gamlss()
bbb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt), 
             method="REML")+pb(income), data=columb)
AIC(b,bb,bbb)
## ploting the fitted models
op<-par(mfrow=c(2,2))
plot(b,scheme=c(0,1))
plot(getSmo(bb), scheme=c(0,1))
par(op)
plot(getSmo(bbb, which=2))
## plot fitted values by district
op<- par(mfrow=c(1,2))
fv <- fitted(b)
names(fv) <- as.character(columb$district)
fv1 <- fitted(bbb)
names(fv1) <- as.character(columb$district)
polys.plot(columb.polys,fv)
polys.plot(columb.polys,fv1)
par(op)}
## bam 
}
\keyword{regression}
