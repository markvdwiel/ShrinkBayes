\name{SpikeSlabPrior}
\alias{SpikeSlabPrior}
\title{
Estimates the fraction of null-hypotheses using a spike-and-slab prior
}
\description{
Estimates the fraction of null-hypotheses using a spike-and-slab prior. Conservative option for inference. 
}
\usage{
SpikeSlabPrior(fitall,fitall0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}.
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

}
\details{
Explicitly maximizes the total marginal likelihood as a function of the fraction of null-hypotheses. 
}
\value{
Numeric
}


\author{
Mark A. van de Wiel
}


\seealso{\code{\link{SpikeSlabUpdatePosterior}} for computing posteriors.
}

