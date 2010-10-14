\name{nn}
\Rdversion{1.1}
\alias{nn}
\alias{nn.control}

\title{A interface function to use  nnet() function within GAMLSS 
}
\description{
The \code{nn()} function is a additive function to be used for GAMLSS models. 
It is an interface for the \code{nnet()} function of package 
\code{nnet} of Brian Ripley. The function \code{nn()} allows the user to use neural networks 
within \code{gamlss}. The great advantage of course comes from the fact GAMLSS models provide a variety of distributions and diagnostics.   
}
\usage{
nn(formula, control = nn.control(...), ...)
nn.control(size = 3, linout = TRUE, entropy = FALSE, softmax = FALSE, 
           censored = FALSE, skip = FALSE, rang = 0.7, decay = 0, 
           maxit = 100, Hess = FALSE, trace = FALSE, 
           MaxNWts = 1000, abstol = 1e-04, reltol = 1e-08)
}
\arguments{
  \item{formula}{ A formula containing the expolanatory variables i.e. ~x1+x2+x3.}
  \item{control}{control to pass the arguments for the nnet() function}
  \item{\dots}{for extra arguments}
  \item{size}{number of units in the hidden layer. Can be zero if there are skip-layer units}
  \item{linout}{switch for linear output units. Default is TRUE, identily link}
  \item{entropy}{switch for entropy (= maximum conditional likelihood) fitting. Default by least-squares.}
  \item{softmax}{switch for softmax (log-linear model) and maximum conditional likelihood fitting. linout, entropy, softmax and censored are mutually exclusive.}
  \item{censored}{A variant on softmax, in which non-zero targets mean possible classes. Thus for softmax a row of (0, 1, 1) means one example each of classes 2 and 3, but for censored it means one example whose class is only known to be 2 or 3.}
  \item{skip}{switch to add skip-layer connections from input to output}
  \item{rang}{Initial random weights on \code{[-rang, rang]}. Value about 0.5 unless the inputs are large, in which case it should be chosen so that \code{rang * max(|x|)} is about 1}
  \item{decay}{parameter for weight decay. Default 0. }
  \item{maxit}{parameter for weight decay. Default 0. }
  \item{Hess}{If true, the Hessian of the measure of fit at the best set of weights found is returned as component Hessian. }
  \item{trace}{switch for tracing optimization. Default FALSE}
  \item{MaxNWts}{The maximum allowable number of weights. There is no intrinsic limit in the code, but increasing MaxNWts will probably allow fits that are very slow and time-consuming. }
  \item{abstol}{Stop if the fit criterion falls below abstol, indicating an essentially perfect fit. }
  \item{reltol}{Stop if the optimizer is unable to reduce the fit criterion by a factor of at least 1 - reltol.}
}

\details{Note that, neural networks are over parameterized models and therefor notorious for multiple maximum.  
There is no guarantee that two identical fits will produce identical results.}
\value{
Note that \code{nn} itself does no smoothing; it simply sets things up for the function \code{gamlss()} which in turn uses the function 
  \code{additive.fit()} for backfitting which in turn uses \code{gamlss.nn()}}
\references{
Rigby, R. A. and  Stasinopoulos D. M. (2005). Generalized additive models for location, scale and shape,(with discussion), 
\emph{Appl. Statist.}, \bold{54}, part 3, pp 507-554.

Ripley, B. D. (1996) Pattern Recognition and Neural Networks. Cambridge. 

Stasinopoulos D. M., Rigby R.A. and Akantziliotou C. (2006) Instructions on how to use the GAMLSS package in R.
Accompanying documentation in the current GAMLSS  help files, (see also  \url{http://www.gamlss.org/}).

Stasinopoulos D. M. Rigby R.A. (2007) Generalized additive models for location scale and shape (GAMLSS) in R.
\emph{Journal of Statistical Software}, Vol. \bold{23}, Issue 7, Dec 2007, \url{http://www.jstatsoft.org/v23/i07}.

Venables, W. N. and Ripley, B. D. (2002) Modern Applied Statistics with S. Fourth edition. Springer. 

}
\author{
Mikis Stasinopoulos \email{d.stasinopoulos@londonmet.ac.uk}, Bob Rigby \email{r.rigby@londonmet.ac.uk}
based on work of  Venables & Ripley wich also based on work by Kurt Hornik and Albrecht Gebhardt.
}


\section{Warning }{You may have to fit the model several time to unsure that you obtain a reasonable minimum} 


\examples{
library(nnet)
data(rock)
area1<- with(rock,area/10000)
peri1<- with (rock,peri/10000) 
rock1<- with(rock, data.frame(perm, area=area1, peri=peri1, shape))
# fit nnet
r1 <- nnet(log(perm)~area+peri+shape, rock1, size=3, decay=1e-3, linout=TRUE, skip=TRUE, max=1000, Hess=TRUE)
summary(r1) 
# get gamlss
library(gamlss) 
cc <- nn.control(size=3, decay=1e-3, linout=TRUE, skip=TRUE, max=1000, Hess=TRUE)
g1 <- gamlss(log(perm)~nn(~area+peri+shape,size=3, control=cc), data=rock1)
summary(g1$mu.coefSmo[[1]])
# predict
Xp <- expand.grid(area=seq(0.1,1.2,0.05), peri=seq(0,0.5, 0.02), shape=0.2)
rocknew <- cbind(Xp, fit=predict(r1, newdata=Xp))
library(lattice)
wf1<-wireframe(fit~area+peri, rocknew, screen=list(z=160, x=-60), aspect=c(1, 0.5), drape=TRUE,  main="nnet()")
rocknew1 <- cbind(Xp, fit=predict(g1, newdata=Xp))
wf2<-wireframe(fit~area+peri, rocknew1, screen=list(z=160, x=-60), aspect=c(1, 0.5), drape=TRUE,  main="nn()")
print(wf1, split=c(1,1,2,1), more=TRUE)
print(wf2, split=c(2,1,2,1))
#--------------------------------------------------------------------------------------- 
 data(rent)
 mr1 <- gamlss(R~nn(~Fl+A, size=5, decay=0.001), data=rent, family=GA)  
 library(gamlss.add)
 mg1<-gamlss(R~ga(~s(Fl,A)), data=rent, family=GA) 
 AIC(mr1,mg1)
newrent <- newrent1 <-newrent2 <- data.frame(expand.grid(Fl=seq(30,120,5), A=seq(1890,1990,5 )))
newrent1$fit <- predict(mr1, newdata=newrent, type="response") ##nn
newrent2$fit <- predict(mg1, newdata=newrent, type="response")# gam
 library(lattice)
 wf1<-wireframe(fit~Fl+A, newrent1, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="nn()")
 wf2<-wireframe(fit~Fl+A, newrent2, aspect=c(1,0.5), drape=TRUE, colorkey=(list(space="right", height=0.6)), main="ga()")
print(wf1, split=c(1,1,2,1), more=TRUE)
print(wf2, split=c(2,1,2,1))
#----------------------------------------------------------------------------------------
data(db)
mdb1 <- gamlss(head~nn(~age,size=20, decay=0.001), data=db)
plot(head~age, data=db)
points(fitted(mdb1)~db$age, col="red")

# do not run
#mdb2 <- gamlss(head~nn(~age,size=20, decay=0.001), data=db, family=BCT)
#plot(head~age, data=db)
#points(fitted(mdb2)~db$age, col="red")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{regression}
