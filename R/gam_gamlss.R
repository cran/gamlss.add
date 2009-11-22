# library(mgcv)
# fit smoothing terms using the gam function of mgcv 
# which is used in the backfitting 
#----------------------------------------------------------------------------------------
ga <-function(formula, ...) 
{ 
#------------------------------------------
# function starts here
#------------------------------------------
    scall <- deparse(sys.call())
if (!is(formula, "formula")) stop("formula argument in ga() needs a formula starting with ~")
# get where "gamlss" is in system call
# it can be in gamlss() or predict.gamlss()  
    rexpr <- grepl("gamlss",sys.calls()) ## 
for (i in 1:length(rexpr))
   { 
 position <- i # get the position
 if (rexpr[i]==TRUE) break
   }
gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
##---
## get the data
if (sys.call(position)[1]=="predict.gamlss()")
     {
      Data <- get("data", envir=gamlss.env)
     }
else { 
      if (is.null(get("gamlsscall", envir=gamlss.env)$data)) stop("the option data in gamlss() is required for ga() to work")
     Data <- get("gamlsscall", envir=gamlss.env)$data
     }
     Data <- data.frame(eval(substitute(Data))) 
      len <- dim(Data)[1] # get the lenth of the data
## out
     xvar <- rep(0,  len) # model.matrix(as.formula(paste(paste("~",ff, sep=""), "-1", sep="")), data=Data) #
   attr(xvar,"formula")     <- formula
   attr(xvar, "gamlss.env") <- gamlss.env
   attr(xvar, "data")       <- as.data.frame(Data)
   attr(xvar, "call")       <- substitute(gamlss.ga(data[[scall]], z, w, ...)) 
   attr(xvar, "class")      <- "smooth"
   xvar
}
#----------------------------------------------------------------------------------------
# the definition of the backfitting additive function
gamlss.ga <-function(x, y, w, xeval = NULL, ...)
{              
   formula <- attr(x,"formula")
   formula <- as.formula(paste("y",deparse(formula), sep=""))
#gamlss.env <- as.environment(attr(x, "gamlss.env"))
      OData <- attr(x,"data") 
      Data <-  if (is.null(xeval)) OData #the trick is for prediction
               else  OData[seq(1,length(y)),]
      Data <- data.frame(eval(substitute(Data)),y,w)
       fit <- gam(formula, data=Data, weights=w, ...)
        df <- sum(fit$edf)-1 
        fv <- fitted(fit) 
 residuals <- y-fv
  if (is.null(xeval))
    {
   list(fitted.values=fv, residuals=residuals,
     nl.df = df, lambda=fit$sp[1], ## we nead df's here 
     coefSmo = fit, var=NA)    # var=fv has to fixed
    }
else 
    {
   ll<-dim(OData)[1]
   pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
    }         

}
      
