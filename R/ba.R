# This is the interface to use the gam() and bam() function of Simon Wood 
# where fitting and initialization are done separately
# Authors: Daniil Kiose and Vlasiis Voudouris
# latest 03-OCT-16 MS
#--------------------------------------------------------------------------------------
ba <-function(formula, control = ba.control(...), ...)  
  {
  #------------------------------------------
  # function starts here
  #------------------------------------------    
  scall <- Reduce(paste, deparse(sys.call(), width.cutoff = 500L)) 
  if (!is(formula, "formula")) 
    stop("formula argument in ba() needs a formula starting with ~")
  # get where "gamlss" is in system call, it can be in gamlss() or predict.gamlss()
  rexpr <- grepl("gamlss",sys.calls()) ## 
  for (i in length(rexpr):1) { 
    position <- i # get the position
    if (rexpr[i]==TRUE) break
  }
  gamlss.env <- sys.frame(position) #gamlss or predict.gamlss
  ## get the data
  ## this has been modified on the 12-12-14 to make sure that 
  ##  if model.frame.gamlss() is used as for example in term.plot() the
  ## function does not fail (It need to implemented to all smoother using formulea?)
  if (sys.call(position)[1]=="predict.gamlss()") { # if predict is used 
    Data <- get("data", envir=gamlss.env)
  } else if (sys.call(position)[1]=="gamlss()") { # if gamlss() is used
    if (is.null(get("gamlsscall", envir=gamlss.env)$data)) { # if no data argument but the formula can be interpreted
      Data <- model.frame(formula)	
    } else {# data argument in gamlss 
      Data <- get("gamlsscall", envir=gamlss.env)$data
    }
  } else {
    Data <- get("data", envir=gamlss.env)
  }
  Data <- data.frame(eval(substitute(Data)))
#-------------------------------------------------
  # new Daniil and Vlasis
  # Initialize bam
  formula <- as.formula(paste0("Y.var", Reduce(paste,deparse(formula, width.cutoff = 500L)) ))
  Data$Y.var = rep(0, nrow(Data))
  #browser()
    G = bam(formula, 
          data = Data,
        offset = control$offset,
        method = control$method,
       control = control$control,
        select = control$select, 
         scale = control$scale, 
         gamma = control$gamma, 
         knots = control$knots, 
            sp = control$sp,
        min.sp = control$min.sp,
       paraPen = control$paraPen, 
    chunk.size = control$chunk.size, 
           rho = control$rho, 
      AR.start = control$AR.start, 
      discrete = control$discrete,
       cluster = control$cluster, 
      nthreads = control$nthreads,
      gc.level = control$gc.level, 
      use.chol = control$use.chol,
       samfrac = control$samfrac, 
          coef = control$coef, 
  drop.unused.levels = control$drop.unused.levels, 
drop.intercept = control$drop.intercept,
             G = NULL,
           fit = FALSE
  )
##bam(formula, family=gaussian(), 
##      data=list()#, weights=NULL, subset=NULL,
#    na.action=na.omit, offset=NULL#, method="fREML"#,control=list()#,
#    select=FALSE#, scale=0#,gamma=1,knots=NULL,sp=NULL,min.sp=NULL,
#    paraPen=NULL,chunk.size=10000,rho=0,AR.start=NULL,discrete=FALSE,
#    cluster=NULL,nthreads=1,gc.level=1,use.chol=FALSE,samfrac=1,
#    coef=NULL,drop.unused.levels=TRUE,G=NULL,fit=TRUE,drop.intercept=NULL,...)
#
 #-------------------------------------------------- 
  xvar <- rep(0,  dim(Data)[1]) 
  attr(xvar,"formula")     <- formula
  attr(xvar,"control")     <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  attr(xvar, "data")       <- as.data.frame(Data)
  attr(xvar, "call")       <- substitute(gamlss.ba(data[[scall]], z, w, ...)) 
  attr(xvar, "class")      <- "smooth"
  attr(xvar, "G")          <- G
  xvar
}

#--------------------------------------------------------------------------------------
ba.control = function(offset = NULL, 
                      method = "fREML", 
                     control = list(), 
                      select = FALSE,
                       scale = 0,
                       gamma = 1, 
                       knots = NULL,
                          sp = NULL, 
                      min.sp = NULL, 
                     paraPen = NULL,
                  chunk.size = 10000,
                         rho = 0,
                    AR.start = NULL,
                    discrete = TRUE,
                     cluster = NULL,
                    nthreads = 2,
                    gc.level = 1,
                    use.chol = FALSE,
                     samfrac = 1,
                        coef = NULL,
          drop.unused.levels = TRUE,
              drop.intercept = NULL, 
                               ...)  
{
  #gam()
  control <- gam.control(...)
  #ga()
  list( offset=offset, method=method, control=control,  select=select,
        scale=scale, gamma=gamma, knots=knots, sp=sp, min.sp=min.sp,
        paraPen = paraPen, chunk.size = chunk.size, rho = rho,
        AR.start = AR.start,
        discrete=discrete,  cluster=cluster, nthreads=nthreads,
        gc.level=gc.level, use.chol=use.chol, samfrac=samfrac,
        coef= coef,
        drop.unused.levels = drop.unused.levels,
        drop.intercept=drop.intercept, ...)
}
#--------------------------------------------------------------------
#--------------------------------------------------------------------------------------
gamlss.ba <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
           Y.var <- y
           W.var <- w
             G <- attr(x,"G")
             control = attr(x,"control")
             G$y <- Y.var
             G$w <- W.var
      G$mf$Y.var <- Y.var
G$mf$`(weights)` <- W.var
             fit <-  bam(G=G, fit=TRUE,
                     offset=control$offset, method=control$method,
                     control=control$control, select=control$select, 
                     scale=control$scale, gamma=control$gamma, 
                     knots=control$knots, sp=control$sp, min.sp=control$min.sp,
                     paraPen=control$paraPen, chunk.size=control$chunk.size, 
                     rho=control$rho, AR.start=control$AR.start, 
                     discrete=control$discrete,
                     cluster=control$cluster, nthreads=control$nthreads,
                     gc.level=control$gc.level, use.chol=control$use.chol,
                     samfrac=control$samfrac, 
                     drop.unused.levels=control$bam$drop.unused.levels)
              df <- sum(fit$edf)-1 
              fv <- fitted(fit) 
       residuals <- y-fv
    list( fitted.values=fv, residuals=residuals,
         nl.df = df, lambda=fit$sp[1], #
         coefSmo = fit, var=NA)    # var=fv has to fixed
  } else { # predict 
      gamlss.env <- as.environment(attr(x, "gamlss.env"))
             obj <- get("object", envir=gamlss.env ) # get the object from predict
              TT <- get("TT", envir=gamlss.env ) # get wich position is now
              SL <- get("smooth.labels", envir=gamlss.env) # all the labels of the smoother
             fit <- eval(parse(text=paste("obj$", get("what", envir=gamlss.env), ".coefSmo[[",as.character(match(TT,SL)), "]]", sep="")))
           OData <- attr(x,"data") 
              ll <- dim(OData)[1]
            pred <- predict(fit,newdata = OData[seq(length(y)+1,ll),])
  }         
}





