\name{NonParaUpdatePosterior}
\alias{NonParaUpdatePosterior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Updates posteriors from a given nonparametric prior, which may include a point mass
}
\description{
This function re-computes posteriors for (one class of) parameters from posteriors obtained by \code{\link{FitAllShrink}} and a new nonparametric prior
obtained by \code{\link{NonParaUpdatePrior}}. 
}
\usage{
NonParaUpdatePosterior(fitall, updateoutput, fitall0=NULL, ncpus = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}.
}
 \item{updateoutput}{
A list object resulting from \code{\link[ShrinkBayes]{NonParaUpdatePrior}}.
}
\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

\item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}
}
\details{
Rescaling of posteriors is used as described in Van de Wiel et al. (2012). 
}
\value{
A list object with the same number of components as the first component of \code{fitall} (number of fits), 
each containing 3-component lists which contain
  \item{postbetanon0}{List of posteriors (one for each parameter/contrast involved)}
  \item{postbeta0}{Point mass (often zero) mixture proportion}
  \item{loglik}{Marginal log-likelihood}
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}

\author{
Mark A. van de Wiel
}

\note{
The resulting posteriors are for the main parameter or contrasts of interest only, which should be indicated in the 
\code{shrinkpara} and \code{shrinklc} option in 
\code{\link{NonParaUpdatePrior}}. Posteriors of other parameters do not alter with respect to those in \code{fitall}.
}

\seealso{\code{\link{NonParaUpdatePrior}} for finding the optimal mixture prior and \code{\link[ShrinkBayes]{FitAllShrink}} for 
fitting under standard parametric priors. In addition, see \code{\link{MixtureUpdatePosterior}} for posteriors given a parametric mixture prior.
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}


