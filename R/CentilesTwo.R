centilesTwo <- function(object, grid.x1, grid.x2, x1.name, x2.name, 
                            cent = 0.05,
                            dist = 0.01, 
                          points = TRUE, 
                           other = list(),
                       point.col = 1, 
                       point.pch = ".", 
                           image = FALSE,
                       image.col = heat.colors(12),       
                          ...)
{
  if (!is.gamlss(object)) stop("the object should be a gamlss object")
  if (is.null(object$call[["data"]])) stop("the gamlss object needs a data argument")
    DaTa <- eval(object$call[["data"]])
     nx1 <- deparse(substitute(x1.name))
     nx2 <- deparse(substitute(x2.name))
  if (!nx1%in%names(DaTa)) stop( "x1.name not in the data")
  if (!nx2%in%names(DaTa)) stop( "x2.name not in the data")
  # get the q function
     fam <- as.gamlss.family(object$call[["family"]])
    npar <- fam$nopar
   fname <- fam$family[[1]]
    qfun <- paste("q", fname, sep = "")
  invcdf <- eval(parse(text=qfun))
  #      IR <- invcdf(.75,...)-invcdf(.25,...
      n1 <- length(grid.x1)
      n2 <- length(grid.x2)
 newdata <- expand.grid(x1=grid.x1, x2=grid.x2)
names(newdata) <- c(nx1, nx2)
  if (length(other))
  {
  newdata <- data.frame(newdata, other[1])
  }
      M3P <- predictAll(object, newdata=newdata, type="response")
       tf <- exclude.too.far(g1=newdata[[nx1]],g2=newdata[[nx2]], d1=DaTa[[nx1]], d2=DaTa[[nx2]],dist=dist)
     PRED <- switch(npar, invcdf(cent, mu=M3P$mu),
                 invcdf(cent, mu=M3P$mu, sigma=M3P$sigma),
                 invcdf(cent, mu=M3P$mu, sigma=M3P$sigma, nu=M3P$nu),
                 invcdf(cent, mu=M3P$mu, sigma=M3P$sigma, nu=M3P$nu, tau=M3P$tau))
  PRED[tf] <- NA
         M <- matrix(PRED, nrow=n1, ncol=n2)
  if (image)
   {
    xlim.image <- ylim.image <- NULL
          eArg <- list(...)
    xlim.image <- ("xlim"%in%names(eArg))*eArg$xlim 
    ylim.image <- ("ylim"%in%names(eArg))*eArg$ylim              
    image(grid.x1, grid.x2, M, col=image.col, ylim=ylim.image, 
          xlim=xlim.image)
    contour(grid.x1, grid.x2, M, add=TRUE, ...) 
   } else   contour(grid.x1, grid.x2, M, ...) 
      
  if (points) points(DaTa[[nx1]], DaTa[[nx2]], pch = point.pch, col=point.col)
}

