\name{ga}
\Rdversion{1.1}
\alias{ga}
\alias{ga.control}

\title{A interface function to use Simon Wood's gam() function within GAMLSS}
\description{
The \code{ga()} function is a additive function to be used within GAMLSS models. 
It is an interface for the \code{gam()} function of package  \code{mgcv} of Simon Wood. The function \code{ga()} allows the user to use all the available smoothers of \code{gam()} within \code{gamlss}. The great advantage of course come from fitting models outside the 
exponential family.    
}
\usage{
ga(formula, control = ga.control(...), ...)

ga.control(offset = NULL, method = "GCV.Cp", optimizer = c("outer", "newton"), 
        control = list(), select = FALSE, knots = NULL,
        sp = NULL, min.sp = NULL, H = NULL, gamma = 1, 
        fit = TRUE, paraPen = NULL, G = NULL, in.out = NULL, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{ A formula containing \code{s()} and \code{te} functions i.e. ~s(x1)+ te(x2,x3).}
  \item{control}{this allow to specify argument within the function \code{gam()} of \code{mgcv}}
  \item{\dots}{arguments used by the \code{gam()} function.}
  \item{offset}{the \code{offset} argument in \code{gam()}}
  \item{method}{the \code{method} argument in \code{gam()}}
  \item{optimizer}{the method optimizer in \code{gam()}}
  \item{select}{the \code{select} argument in \code{gam()}}
  \item{knots}{the \code{knots} argument in \code{gam()}}
  \item{sp}{the \code{sp} argument in \code{gam()}}
  \item{min.sp}{
the \code{min.sp} argument in \code{gam()}
}
  \item{H}{
the \code{H} argument in \code{gam()}
}
  \item{gamma}{
the \code{gamma} argument in \code{gam()}
}
  \item{fit}{
the \code{fit} argument in \code{gam()}
}
  \item{paraPen}{
the \code{paraPen} argument in \code{gam()}
}
  \item{G}{
the \code{G} argument in \code{gam()}
}
  \item{in.out}{
the \code{in.out} argument in \code{gam()}
}

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
## normal error additive in Fl and A
ga2 <- gam(R~s(Fl)+s(A), method="REML", data=rent)
gn2 <- gamlss(R~ga(~s(Fl)+s(A), method="REML"),  data=rent) # additive
gb2 <- gamlss(R~pb(Fl)+pb(A), data=rent) # additive
AIC(ga2,gn2, gb2, k=0)
AIC(ga2,gn2, gb2)## normal error additive in Fl and A
ga2 <- gam(R~s(Fl)+s(A), method="REML", data=rent)
gn2 <- gamlss(R~ga(~s(Fl)+s(A), method="REML"),  data=rent) # additive
gb2 <- gamlss(R~pb(Fl)+pb(A), data=rent) # additive
AIC(ga2,gn2, gb2, k=0)
AIC(ga2,gn2, gb2)
#---------------------------------------------------------
## gamma error two variables s() f unction
ga4 <-gam(R~s(Fl,A), method="REML", data=rent, family=Gamma(log))
gn4 <- gamlss(R~ga(~s(Fl,A), method="REML"), data=rent, family=GA) 
AIC(ga4,gn4, k=0)
AIC(ga4,gn4)
op<-par(mfrow=c(1,2))
vis.gam(ga4)
vis.gam(gn4$mu.coefSmo[[1]])
par(op)
#---------------------------------------------------------
## predict
newrent <- data.frame(expand.grid(Fl=seq(30,120,5), A=seq(1890,1990,5 )))
newrent1 <-newrent2 <- newrent
newrent1$pred <- predict(ga4, newdata=newrent, type="response")
newrent2$pred <- predict(gn4, newdata=newrent, type="response")
library(lattice)
wf1<-wireframe(pred~Fl*A, newrent1, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gam()")
wf2<-wireframe(pred~Fl*A, newrent2, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="gamlss()")
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
vis.gam(gn5$mu.coefSmo[[1]])
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
bb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt), method="REML"), data=columb)
AIC(b,bb, k=0)
plot(b,scheme=1)
plot(bb$mu.coefSmo[[1]], scheme=1)
par(op)
## An important covariate added...
b <- gam(crime ~ s(district,bs="mrf",k=20,xt=xt)+s(income),
         data=columb,method="REML")
## x in gam() 
bb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt)+s(income), method="REML"), data=columb)
## x in gamlss()
bbb <- gamlss(crime~ ga(~s(district,bs="mrf",k=20,xt=xt), method="REML")+pb(income), data=columb)
AIC(b,bb,bbb)

op<-par(mfrow=c(1,2))
plot(b,scheme=c(0,1))
par(op)

op<-par(mfrow=c(1,2))
plot(bb$mu.coefSmo[[1],scheme=c(0,1))
par(op)

plot(bbb$mu.coefSmo[[2]])

## plot fitted values by district
op<- par(mfrow=c(1,2))
fv <- fitted(b)
names(fv) <- as.character(columb$district)
fv1 <- fitted(bbb)
names(fv1) <- as.character(columb$district)
polys.plot(columb.polys,fv)
polys.plot(columb.polys,fv1)
par(op)}

}
\keyword{regression}
