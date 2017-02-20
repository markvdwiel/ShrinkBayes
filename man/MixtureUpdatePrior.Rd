\name{MixtureUpdatePrior}
\alias{MixtureUpdatePrior}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Replacing a non-mixture prior by an optimal, symmetric mixture prior
}
\description{
This is the direct maximization procedure, outlined in the Van de Wiel et al. (2012) reference below. It allows to replace one non-mixture 
prior by a symmetric mixture prior, which may be desirable for the main parameter of interest. 
}
\usage{
MixtureUpdatePrior(fitall,fitall0=NULL, shrinkpara=NULL, modus="mixt", shrinklc=NULL,  lincombs=NULL,ntotal = 10000, maxsupport=6,pointmass=0,
pminvec = c(0,0.25,0.5,0.75,1),p0vec = c(0.5,0.7,0.9,0.95,0.99,0.999,1), 
meanvec = c(0.1, 0.3, 0.5, 0.75,1.5),sdvec=c(0.2,0.5,0.8,1.5,3),meansdauto=TRUE,
ncpus=2,refinegrid=TRUE,symmetric=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

  \item{shrinkpara}{
Character or character string. Name(s) of the variable(s) for which the mixture prior is fit. Corresponding variable can be a factor. }

\item{modus}{
Character string. Parametric form of the continuous component. Either \code{"mixt"}, 
\code{"gauss"} or \code{"laplace"} for 2-component Gaussian mixture, Gaussian, and Laplace, respectively. 
}

 \item{shrinklc}{
Character string. Name of the linear combination for which the nonparametric prior is fit.
}

\item{lincombs}{
List object. Name of the list object that contains the linear combination(s), usually created by \code{inla.make.lincomb} or \code{AllComp}. 
Only required when \code{fitall} does not consist of INLA outputs, e.g. after \code{CombinePosteriors}.
}
  
  \item{ntotal}{
Integer. Number of posteriors that are used to determine the new prior.
}
  \item{maxsupport}{
Numeric. maximum of the support of the prior and posteriors. For numerical stability. -maxsupport is the minimum of the support. 
}

\item{pointmass}{
Numeric. Location of the pointmass. 
}

\item{pminvec}{
Numerical vector. Grid values for probability mass on the negative component with respect to the continuous component. 
Only relevant for \code{modus="mixt"}.
} 

  \item{p0vec}{
Numerical vector. Grid values for probability mass on the point mass.
}

  \item{meanvec}{
Numerical vector. Grid values for mean of the continous component.
}

  \item{sdvec}{
Numerical vector. Grid values for standard deviatioon of the continous component.
}

\item{meansdauto}{
Boolean. If \code{TRUE}: automatically computes initial grid for \code{meanvec} and \code{sdvec}. Overrules setting 
for those two parameters.
}
  \item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}

\item{refinegrid}{
Boolean. If \code{TRUE}: automatically refines grid for parameters once the initial optimal values are found. 
}

\item{symmetric}{
Boolean. If \code{TRUE}: forces a symmetrix prior. Only relevant for mixture priors 
}
}
\details{
This function corresponds to the direct maximization procedure in Van de Wiel et al. (2012). The procedure is currently only implemented for fixed 
regression parameters or functions thereof. Also, only symmetric priors are currently supported.  
About \code{shrinklc}: it is assumed that only one type of linear combinations is present in the fit object \code{fitall}.
}
\value{
A list with two components
  \item{allparam}{Numerical matrix with rows ordered according to log marginal likelihood, containing parameter values of the mixture prior and log marginal likelihood}
  \item{inputpar}{List with input parameters used}
}

\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}

\author{
Mark A. van de Wiel
}

\note{
Computing time increases proportionally with the product of the length of the parameters \code{p0vec, meanvec, sdvec} and 
of \code{p0widevec, sdwidevec} if \code{addwide=TRUE}. After a first run, it may be good to do a second one on a finer grid. 
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{NonParaUpdatePrior}} for non-parametric priors and \code{\link{MixtureUpdatePosterior}} for computing posteriors from the output
of this function. 
}

