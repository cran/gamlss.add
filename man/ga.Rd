\name{ga}
\Rdversion{1.1}
\alias{ga}
\title{A interface function to use Simon Wood's gam() function within GAMLSS }
\description{
The \code{ga()} function is a additive function to be used for GAMLSS models. 
It is an interface for the \code{gam()} function of package 
\code{mgcv} of Simon Wood. The function \code{ga()} allows the user to use all of the available smoothers 
of \code{gam()} within \code{gamlss}. The great advantage of course come from fitting models outside the 
exponential family.    
}
\usage{
ga(formula, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ A formula containing \code{s()} and \code{te} functions i.e. ~s(x1)+ te(x2,x3).}
  \item{\dots}{arguments used by the \code{gam()} function.}
}
\details{
Note that \code{ga} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function 
  \code{additive.fit()} for backfitting which in turn uses \code{gamlss.ga()}

Note that, in our (limited) experience,  for normal errors or exponential family, the  fitted models using \code{gam()} 
and \code{ga()} within \code{gamlss()} are identical
or at least very similar.  This is particularly true if the default values for \code{gam()} are used.   
   
}
\value{
the fitted values of the smoother is returned, endowed with a number of attributes. 
The smoother fitted values are used in the construction of the overall fitted values of the particular distribution parameter.
The attributes can be use to obtain information about the individual fit. In particular the \code{coefSmo} within the parameters
of the fitted model contains the final additive fit.
}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Wood S.N. (2006) Generalized Additive Models: An Introduction with R. Chapman and Hall/CRC Press. 
}
\author{Mikis Stasinopoulos}

\section{Warning }{The function id experimental so please report any peculiar behaviour to the authors} 


\examples{
library(gamlss)
library(mgcv)
data(rent)
#---------------------------------------------------------
# normal errors one x-variable
ga1 <-  gam(R~s(Fl, bs="ps", k=20), data=rent, method="ML")
gn1 <- gamlss(R~pb(Fl), data=rent) # additive
gn2 <- gamlss(R~ga(~s(Fl)), data=rent) # additive
AIC(ga1,gn1,gn2,  k=0)
#--------------------------------------------------------
# normal error additive in Fl and A
# normal error additive in Fl and A
ga2 <- gam(R~s(Fl)+s(A), data=rent)
gn0 <- gamlss(R~pb(Fl)+pb(A), data=rent) # additive
gn1 <- gamlss(R~ga(~s(Fl)+s(A)), data=rent) # additive
gn2 <- gamlss(R~ga(~s(Fl))+ga(~s(A)), data=rent)
# similar fitting 
AIC(ga2,gn0,gn1, gn2,  k=0)
#---------------------------------------------------------
# gamma errors one x-var
ga1 <-    gam(R~s(Fl), data=rent, family=Gamma)
gg1 <- gamlss(R~pb(Fl), data=rent, family=GA)
gg2 <- gamlss(R~ga(~s(Fl)), data=rent, family=GA)
AIC(ga1, gg1, gg2,  k=0)
# different degrees of freedon for mu in ga1
#---------------------------------------------------------
# gamma error two variables s() function
g22 <-gam(R~s(Fl,A), data=rent, family=Gamma)
gm22 <- gamlss(R~ga(~s(Fl,A)), data=rent, family=GA) 
AIC(g22,gm22)
# predict
newrent <- data.frame(expand.grid(Fl=seq(30,120,5), A=seq(1890,1990,5 )))
newrent$pred2 <- predict(gm22, newdata=newrent, type="response")
newrent$pred1 <- predict(g22, newdata=newrent, type="response")
library(lattice)
wf1<-wireframe(pred1~Fl*A, newrent, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gam()")
wf2<-wireframe(pred2~Fl*A, newrent, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gamlss()")
print(wf1, split=c(1,1,2,1), more=TRUE)
print(wf2, split=c(2,1,2,1))
#---------------------------------------------------------
#gamma error two variables te() function
g221 <-gam(R~te(Fl,A), data=rent, family=Gamma)
gm221 <- gamlss(R~ga(~te(Fl,A)), data=rent, family=GA) 
AIC(g221,gm221)
# predict
newrent <- data.frame(expand.grid(Fl=seq(30,120,5), A=seq(1890,1990,5 )))
newrent$pred2 <- predict(gm221, newdata=newrent, type="response")
newrent$pred1 <- predict(g221, newdata=newrent, type="response")
library(lattice)
wf1<-wireframe(pred1~Fl*A, newrent, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gam()")
wf2<-wireframe(pred2~Fl*A, newrent, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gamlss()")
print(wf1, split=c(1,1,2,1), more=TRUE)
print(wf2, split=c(2,1,2,1))
#----------------------------------------------------------

}
\keyword{regression}
