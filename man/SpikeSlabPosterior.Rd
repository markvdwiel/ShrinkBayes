\name{SpikeSlabUpdatePosterior}
\alias{SpikeSlabUpdatePosterior}
\title{
Updates posteriors from a fit-object and a spike-and-slab prior}

\description{
Updates posteriors using a spike-and-slab prior.
}
\usage{
SpikeSlabUpdatePosterior(fitall,fitall0, p0est,  modus="fixed", shrinkpara=NULL, shrinklc=NULL, ncpus=1,only0=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from 
\code{\link[ShrinkBayes]{CombinePosteriors}}.
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

 \item{p0est}{
Numeric. Estimate for fraction of null-hypotheses, usually obtained from \code{\link[ShrinkBayes]{SpikeSlabPrior}}
}

\item{modus}{
Character string. Either \code{"fixed"}, \code{"random"} or \code{"logdisp"}. Type of variable for which the posterior is desired.
}
  \item{shrinkpara}{
Character or character vector. Name(s) of the variable(s) for which the posterior is desired. Corresponding variable can be a factor. 
}

\item{shrinklc}{
Character string. Name of the linear combination for which the posterior is desired.
}

\item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}

\item{only0}{
Boolean.  If \code{TRUE}, the function only computes the posterior probability of the null-model.
}
}
\details{
This function allows one to use spike-and-slab priors, which usually results in a conservative False Discovery Rate estimate. 
}
\value{
A list object with the same number of components as the first component of \code{fitall} (number of fits),
each containing 3-component lists which contain
  \item{postbetanon0}{Continuous component of the posterior (as matrix)}
  \item{postbeta0}{Point mass (often zero) mixture proportion}
  \item{loglik}{Marginal log-likelihood}
}

%\references{
%Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012).
%Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
%}

\author{
Mark A. van de Wiel
}

\seealso{\code{\link{SpikeSlabPrior}} for computing the proportion of null-hypotheses.
}

