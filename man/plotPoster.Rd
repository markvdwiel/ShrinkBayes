\name{plotPoster}
\alias{plotPoster}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots posterior densities.
}
\description{
Plots posterior densities, possibly including a posterior null point mass.
}
\usage{
plotPoster(id, poster, miny = 10^(-4), xlabel = "", ylabel = "", plotlabels = NULL, cols = NULL, ltys = NULL,lwds =NULL, xran = NULL, 
xscatter = NULL, legendxy = NULL, legendbox = "n", include0 = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{id}{
Integer. Index of feature for which the posterior density is desired
}
  \item{poster}{
A list of lists where each sublist contains posteriors of (possibly) multiple parameters (or linear combinations) from the same data.
A posterior is represented as a matrix with two columns denoting support and posterior density.
}
  \item{miny}{
Numeric. Minimum value of density (in the tails)
}
  \item{xlabel}{
Character. Label for the x-axis.
}
  \item{ylabel}{
Character. Label for the y-axis.
}
  \item{plotlabels}{
Character vector. Label(s) for the graphs.
}
  \item{cols}{
Numerical or character vector. Colours of the graphs.
}
  \item{ltys}{
Numerical or character vector. Line types for the graphs.
}

\item{lwds}{
Numerical vector. Line widths for the graphs.
}

  \item{xran}{
Numerical 2-component vector. Range for the x-values.
}
  \item{xscatter}{
Numeric. Scatter for the null point masses.
}
  \item{legendxy}{
Numerical 2-component vector. Determines (x,y)-location of the legend.
}
  \item{legendbox}{
Character. If \code{legendbox="n"} no legend box is drawn; if \code{legendbox="o"} a legend box is drawn.
}
  \item{include0}{
Boolean. \code{include0=TRUE} the null point mass is drawn as well.
}
}
\details{
About \code{posteriorlist}: this usually results from \code{\link{NonParaUpdatePosterior}}, \code{\link{MixtureUpdatePosterior}}, 
\code{\link{BFUpdatePosterior}} or \code{\link{SpikeSlabUpdatePosterior}}. About \code{xscatter}: when multiple densities are displayed 
in one figure, the null point masses overlap. Hence, a small scatter is added to avoid the overlap. 
}
\value{
Plot.
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012).
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}
\author{
Mark A. van de Wiel
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SummaryTable}} for a summary table of all results.
}

