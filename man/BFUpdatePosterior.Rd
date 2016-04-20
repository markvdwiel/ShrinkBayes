\name{BFUpdatePosterior}
\alias{BFUpdatePosterior}
\title{
Updates posteriors for testing nested models that differ by more than one variable or simply extracts all posteriors for one 
particular parameter 
}

\description{
Updates posteriors from output by  \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}. 
Useful for testing nested models or for simply extracting posteriors for a particular parameter after applying 
\code{\link[ShrinkBayes]{FitAllShrink}}.}

\usage{
BFUpdatePosterior(fitall, updateoutput,fitall0=NULL, modus="fixed", shrinkpara=NULL, shrinklc=NULL, ncpus=1,only0=FALSE)
}

\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}.
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

 \item{updateoutput}{
A list object resulting from \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}.
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
This function is mostly useful for testing nested models that differ by more than one variable. 
}
\value{
A list object with the same number of components as the first component of \code{fitall} (number of fits), 
each containing 3-component lists which contain
  \item{postbetanon0}{Continuous component of the posterior (as matrix)}
  \item{postbeta0}{Point mass (often zero) mixture proportion. Zero if \code{fitall0=NULL}}
  \item{loglik}{Marginal log-likelihood}
}

%\references{
%Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
%Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
%}

\author{
Mark A. van de Wiel
}


