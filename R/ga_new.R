# This is the interface to use the gam() and bam() function of Simon Wood 
# where fitting and initialization are done separately
# Authors: Daniil Kiose and Vlasiis Voudouris
# latest 03-Oct-16 MS
#--------------------------------------------------------------------------------------
ga <-function(formula, control = ga.control(...), ...)  
  {
  #------------------------------------------
  # function starts here
  #------------------------------------------    
  scall <- deparse(sys.call(), width.cutoff = 500L) 
  if (!is(formula, "formula")) 
    stop("formula argument in ga() needs a formula starting with ~")
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
  # Initialize gam
     formula <- as.formula(paste0("Y.var", deparse(formula, width.cutoff = 500L)))
  Data$Y.var <- rep(0, nrow(Data))
          G <- gam(formula, 
            data = Data,
          offset = control$offset, 
          method = control$method, 
       optimizer = control$optimizer, 
         control = control$control,
           scale =  control$scale,
          select = control$select, 
           knots = control$knots, 
              sp = control$sp,
          min.sp = control$min.sp, 
               H = control$H, 
           gamma = control$gamma,  
         paraPen = control$paraPen,
          in.out = control$in.out, 
drop.unused.levels = control$drop.unused.levels,
  drop.intercept = control$drop.intercept,
        discrete = control$discrete,
              G = NULL, 
            fit = FALSE)
 #-------------------------------------------------- 
  xvar <- rep(0,  dim(Data)[1]) 
  attr(xvar,"formula")     <- formula
  attr(xvar,"control")     <- control
  attr(xvar, "gamlss.env") <- gamlss.env
  attr(xvar, "data")       <- as.data.frame(Data)
  attr(xvar, "call")       <- substitute(gamlss.ga(data[[scall]], z, w, ...)) 
  attr(xvar, "class")      <- "smooth"
  attr(xvar, "G")          <- G
  xvar
}

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
ga.control = function(offset = NULL, 
                      method = "REML", 
                   optimizer = c("outer","newton"), 
                     control = list(),
                       scale = 0,
                      select = FALSE, 
                       knots = NULL,
                          sp = NULL, 
                      min.sp = NULL, 
                           H = NULL, 
                       gamma = 1, 
                     paraPen = NULL,  
                      in.out = NULL,
          drop.unused.levels = TRUE,
              drop.intercept = NULL,
                    discrete = FALSE,
          ...)
{
  #gam()
  control <- gam.control(...)
  #ga()
  list(offset=offset, method=method, optimizer=optimizer, control=control, 
       scale= scale, 
       select=select, knots=knots, sp=sp, min.sp=min.sp, H=H, gamma=gamma, 
       paraPen=paraPen, in.out=in.out,   drop.unused.levels = drop.unused.levels,
       drop.intercept=drop.intercept, discrete = discrete, ...)
}
#--------------------------------------------------------------------------------------
gamlss.ga <-function(x, y, w, xeval = NULL, ...) {
  if (is.null(xeval))
  {#fitting
   #formula <- attr(x,"formula")
    #control <- as.list(attr(x, "control"))
           Y.var <- y
           W.var <- w
               G <- attr(x,"G") 
             G$y <- Y.var
             G$w <- W.var
      G$mf$Y.var <- Y.var
G$mf$`(weights)` <- W.var
             fit <-  gam(G=G, fit=TRUE)
              df <- sum(fit$edf)-1 
              fv <- fitted(fit) 
       residuals <- y-fv
    list(fitted.values=fv, residuals=residuals,
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





