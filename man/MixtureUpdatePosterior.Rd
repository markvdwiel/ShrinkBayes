\name{MixtureUpdatePosterior}
\alias{MixtureUpdatePosterior}
\title{
Updates posteriors from a given mixture prior
}
\description{
This function re-computes posteriors for (one class of) parameters from posteriors obtained by \code{\link{FitAllShrink}} and a new mixture prior
obtained by \code{\link{MixtureUpdatePrior}}. 
}
\usage{
MixtureUpdatePosterior(fitall, updateoutput, fitall0=NULL, ncpus = 1, robustlarge=TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}.
}

 \item{updateoutput}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{MixtureUpdatePrior}}.
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

\item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}

\item{robustlarge}{
Boolean. Add a small proportion of flat Gaussian prior to robustify results for large effects.
}

}
\details{
Rescaling of posteriors is used as described in Van de Wiel et al. (2012). About \code{robustlarge}: in rare cases 
the fitted prior may mismatch extreme effects. For that we weigh the continuous part of the prior with a very small proportion (0.001) of 
a fairly flat prior (central gaussian with sd=5). This has almost no effect for most posteriors, but renders a better posterior for extreme effects. 
}
\value{
A list object with the same number of components as the first component of \code{fitall} (number of fits), 
each containing 3-component lists which contain
  \item{postbetanon0}{Continuous component of the posterior (as matrix)}
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
\code{\link{MixtureUpdatePrior}}. Posteriors of other parameters do not alter with respect to those in \code{fitall}.
}

\seealso{\code{\link{MixtureUpdatePrior}} for finding the optimal mixture prior and \code{\link[ShrinkBayes]{FitAllShrink}} for 
fitting under standard parametric priors. In addition, see \code{\link{NonParaUpdatePosterior}} for posteriors given a nonparametric prior.
}

