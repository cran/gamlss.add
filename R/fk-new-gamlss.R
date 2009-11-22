# ms Tuesday, July 7, 2009 
# fit smoothing terms using the  curfit.free.knot() function
# which is used in the backfitting 
#----------------------------------------------------------------------------------------
fk <-function(x, degree=3, start=NULL, ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call())
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()       
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in length(rexpr):1)
   { 
 position <- i # get the position, we are geting the fist from the last
 if (rexpr[i]==TRUE) break
   }
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
#--------
## get a random name to use it in the gamlss() environment
#--------
               sl <-sample(letters, 4)
      fourLetters <- paste(paste(paste(sl[1], sl[2], sep=""), sl[3], sep=""),sl[4], sep="")
  startLambdaName <- paste("start.Lambda",fourLetters, sep=".")
## put the starting values in the gamlss()environment
#--------
   assign(startLambdaName, start, envir=gamlss.env)
      len <- length(x) # get the lenth of the data
## out
     xvar <- x#rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
  # attr(xvar,"formula")     <- formula
   attr(xvar, "x")             <- x
   attr(xvar, "gamlss.env")    <- gamlss.env
   attr(xvar, "degree")        <- degree
   attr(xvar, "NameForLambda") <- startLambdaName
   attr(xvar, "call")          <- substitute(gamlss.fk(data[[scall]], z, w, ...)) 
   attr(xvar, "class")         <- "smooth"
   xvar
}
##---------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.fk <-function(x, y, w, xeval = NULL, ...)
{              
      xvar <-  if (is.null(xeval)) attr(x,"x")
             else  attr(x,"x")[seq(1,length(y))]
         degree <- attr(x,"degree")
     gamlss.env <- as.environment(attr(x, "gamlss.env"))
startLambdaName <- as.character(attr(x, "NameForLambda"))  
         lambda <- get(startLambdaName, envir=gamlss.env) ## geting the starting knots 
            fit <- fitFreeKnots(x=xvar, y=y, w=w, degree=degree, knots = lambda,   ...)
       #cat("knot", knots(fit), "\n")
        assign(startLambdaName, fit$breakPoints, envir=gamlss.env)
  if (is.null(xeval))
   {
   list(fitted.values=fitted(fit), residuals=resid(fit),  nl.df = fit$df-2, lambda=fit$knots[1], ## we nead df's here 
     coefSmo = fit, var=NA)
   }    
else 
   {
      lenN <- length(attr(x,"x"))
     nxval <- attr(x,"x")[(length(y)+1):lenN]
      pred <- predict(fit,newdata=nxval)
      pred  
   }         

}
      
