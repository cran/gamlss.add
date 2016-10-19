\name{plot.nnet}
\alias{plot.nnet}

\title{Plotting fitted neural networks
}
\description{
A function to plot the results of a neural network fit based on the \code{plotnet()} function of the package \pkg{NeuralNetTools}
}
\usage{
## S3 method for class 'nnet'
\method{plot}{nnet}(x, nid = TRUE, all.out = TRUE, all.in = TRUE, bias = TRUE, 
wts.only = FALSE, rel.rsc = 5, circle.cex = 5, node.labs = TRUE, 
var.labs = TRUE, x.lab = NULL, y.lab = NULL, line.stag = NULL, 
struct = NULL, cex.val = 1, alpha.val = 1, circle.col = "lightblue",
pos.col = "black", neg.col = "grey", max.sp = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{A neural network fitted model }
  \item{nid}{logical value indicating if neural interpretation diagram is plotted, default is  code{TRUE}}
\item{all.out}{character string indicating names of response variables for which connections are plotted, default all}
 \item{all.in}{character string indicating names of input variables for which connections are plotted, default all}
  \item{bias}{logical value indicating if bias nodes and connections are plotted, not applicable for networks from mlp function, default \code{TRUE}}
  \item{wts.only}{logical value indicating if connections weights are returned rather than a plot, default \code{FALSE}}
  \item{rel.rsc}{numeric value indicating maximum width of connection lines, default 5}
  \item{circle.cex}{numeric value indicating size of nodes, passed to cex argument, default 5}
  \item{node.labs}{logical value indicating if text labels are plotted, default \code{TRUE}}
  \item{var.labs}{logical value indicating if variable names are plotted next to nodes, default \code{TRUE}}
  \item{x.lab}{character string indicating names for input variables, default from model object}
  \item{y.lab}{character string indicating names for output variables, default from  model object}
  \item{line.stag}{numeric value that specifies distance of connection weights from nodes}
 \item{struct}{numeric value of length three indicating network architecture (no nodes for input, hidden, output), required only if mod.in is a numeric vector}
\item{cex.val}{numeric value indicating size of text labels, default 1}
\item{alpha.val}{numeric value (0-1) indicating transparency of connections, default 1}
  \item{circle.col}{text value indicating colour of nodes default "lighrblue" }
 \item{pos.col}{text value indicating colour of the possitive connections, default "black" }
 \item{neg.col}{text value indicating colour of the negative connections, default "gray"}
  \item{max.sp}{logical value indication whether the space betwwen nodes in each laers is maximised}
  \item{\dots}{for further arguments}
}
\details{
The function \code{plot.nnet()} is (almost) identical to the function \code{plot.nnet()} created by Marcus W. Beck it was first published in the web but  now is part of the \pkg{NeuralNetTools} package in \R{} under the name \code{plotnet()}. Here we modify the function it so it works within the \pkg{gamlss.add} package. This involves of borrowing the functions \code{rescale()}, \code{zero_range()} and \code{alpha()} from package \pkg{scales}. 
}
\value{
The function is producing a plot
}
\references{
Marcus W. Beck (2015). NeuralNetTools: Visualization and Analysis Tools for Neural Networks. R package version 1.4.1. \url{https://cran.r-project.org/package=NeuralNetTools}

Hadley Wickham (2014). scales: Scale functions for graphics. R package version 0.4.0. \url{https://cran.r-project.org/package=scales}

}
\author{Marcus W. Beck <mbafs2012@gmail.com> modified by Mikis Stasinopoulos}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{nn}}}
\examples{
r1 <- gamlss(R~nn(~Fl+A+H+loc, size=10, decay=0.2), data=rent, 
      family=GA, gd.tol=1000, n.cyc=5)
getSmo(r1)
plot(getSmo(r1), y.lab=expression(eta[1]))
plot(getSmo(r1), y.lab=expression(g[1](mu)))
\dontrun{
r2 <- gamlss(R~nn(~Fl+A+H+loc, size=10, decay=0.2), 
      sigma.fo=~nn(~Fl+A+H+loc, size=10, decay=0.2),data=rent, 
      family=GA, gd.tol=1000, n.cyc=5)
plot(getSmo(r2), y.lab=expression(g[1](mu)))
plot(getSmo(r2, what="sigma"), y.lab=expression(g[2](sigma)))
}
}
\keyword{regression}
