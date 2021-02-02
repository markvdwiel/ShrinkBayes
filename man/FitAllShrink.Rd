\name{FitAllShrink}
\alias{FitAllShrink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Applies INLA to multiple data rows using shrinkage priors.
}
\description{
This is a wrapper function that facilitates multiple calls to INLA using parallel processing. It currently supports the following 
univariate likelihoods: Poisson, negative binomial, zero-inflated negative binomial, Gaussian. It uses the output from 
\code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}} as input for the prior parameters. 


}
\usage{
FitAllShrink(forms, dat, shrinksimul, finalprior=FALSE, dispersefixed = 10, disperseaddfixed = 1, disperserandom = 1, maxprecfixed = 4, fams = "zinb", ncpus = 2, 
effoutput = TRUE, keepmargrand = FALSE, keepmarghyper = TRUE, setthreads1 = TRUE, showupdate = FALSE, silentINLA = TRUE, updateby = 5000, 
ndigits = 5, addpackage = NULL, safemode = TRUE, designlist=NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{forms}{
Formula, or list of formulas the length of which equals the number of data rows. See \code{\link[INLA]{inla}} and 
\code{\link[INLA]{f}} for specification of the model formula.
}
  \item{dat}{
Matrix, data frame or list containing the data. Rows are features, columns are samples. For lists: 
each component represents a feature.
}
  \item{shrinksimul}{
A list object resulting from \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}
}

\item{finalprior}{
Boolean. If TRUE, \code{dispersefixed}, \code{disperseaddfixed}, \code{disperserandom} are all set to 1.
}
  \item{dispersefixed}{
Numeric. Inflation factor for the variance of the main fixed parameter.  
}

  \item{disperseaddfixed}{
Numeric. Inflation factor for the variance of the additional fixed parameter.  
}

  \item{disperserandom}{
Numeric. Inflation factor for the variance of the random effects precision.  
}
  \item{maxprecfixed}{
Numeric. Maximum precision used for the main fixed effect.
}
  \item{fams}{Character string. Either equal to \code{"poisson"}, \code{"zip"} 
(zero-inflated Poisson), \code{"nb"} (negative binomial), or \code{"zinb"}
(zero-inflated negative binomial), \code{"gaussian"}: likelihood to be used.

}
  \item{ncpus}{Integer.  The number of cpus to use for parallel computations.
}
  \item{effoutput}{
Boolean. If FALSE, all INLA output will be saved. If TRUE, some fields will be deleted.
}
  \item{keepmargrand}{
Boolean. Do you wish to save the marginals of the random effect regression parameters (beta's)?
}
  \item{keepmarghyper}{
Boolean. Do you wish to save the marginals of the hyper-parameters?
}
  \item{setthreads1}{
Boolean. If TRUE, sequential computation is forced within each call to \code{\link[INLA]{inla}}. 
}
  \item{showupdate}{
Boolean. Do you wish to see updates on progression of the computation? TRUE may slow down the computations due to less efficient 
parallel computation.}
  \item{silentINLA}{
Boolean. Do you wish to silence the output of \code{\link[INLA]{inla}}? NOTE: currently not all INLA output can be silenced. 
}
  \item{updateby}{
Integer, only relevant when showupdate=TRUE. Show an update for each \code{updateby} number of features executed.
}
  \item{ndigits}{
Integer. Numerical precision in digits for the output.
}
  \item{addpackage}{
Character string. Additional package that you wish to export to slave nodes when parallel computing.
}
  \item{safemode}{
Boolean. Only relevant for \code{fams="zinb"} or \code{fams="nb"}. If fitting fails, should it re-run in a safe mode?
}

\item{designlist}{
List. Components are data frames containting the variables in \code{form}. Length of list should equal nr of features in \code{data}. If \code{NULL} design is assumed to be the same for all features (and inferred from \code{form})
}
  
  \item{\dots}{
Further arguments passed on to \code{\link[INLA]{inla}}.
}
}
\details{ 
\code{dispersefixed}, \code{disperseaddfixed} and \code{disperserandom} can be used to fit under a flatter prior than the one found by the joint shrinkage procedures
\code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}. This is typically useful when one aims to empirically fit 
a mixture prior or nonparametric prior for the main fixed or random parameter (see \code{\link[ShrinkBayes]{MixtureUpdatePrior}} or 
\code{\link[ShrinkBayes]{NonParaUpdatePrior}}). \code{finalprior=FALSE} should be used in combination with \code{\link[ShrinkBayes]{MixtureUpdatePrior}} or 
\code{\link[ShrinkBayes]{NonParaUpdatePrior}}, while \code{finalprior=TRUE} should be used in combination with 
\code{\link[ShrinkBayes]{BFUpdatePosterior}} or \code{\link[ShrinkBayes]{SpikeSlabUpdatePosterior}}, because the latter two functions do 
not further update the prior. 

About \code{setthreads1}: usually it is computationally most efficient to use only one node (thread) per INLA call. About 
\code{addpackage}: sometimes the formula contains a call to an external package (e.g. for fitting splines). If so, this package needs to be specified for 
exporting to slave nodes when computing parallelly. About \code{safemode}: if TRUE, then features for which fitting fails when the size parameter of the NB 
or ZI-NB are estimated are re-fitted with fixed size parameter, which equals the prior mean. 

}

\value{
Two-component list object.
  \item{res}{A list of length \code{nrow(dat)} containing \code{\link[INLA]{inla}} output for each data row}
  \item{priors}{A list passing through all information on priors present in \code{shrinksimul}}
}

\references{
Rue H, Martino S and Chopin N (2009). Approximate Bayesian inference for latent
Gaussian models by using integrated nested Laplace approximations (with discussion). J. R.
Statist. Soc. B, 71, 319-392. www.r-inla.org

Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}

\author{
Mark A. van de Wiel
}
\note{
Extensions to the supported likelihoods (including multivariate ones) may be released in the feature
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\url{www.r-inla.org}, \code{\link[ShrinkBayes]{ShrinkSeq}}, \code{\link[ShrinkBayes]{ShrinkGauss}}
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}
