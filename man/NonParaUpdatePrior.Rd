\name{NonParaUpdatePrior}
\alias{NonParaUpdatePrior}
\title{
Replacing a parametric prior by an optimal non-parametric one}
\description{
This is the iterative marginal procedure, outlined in the Van de Wiel et al. (2012) reference below. It allows to replace one parametric 
prior by a nonparametric one, which is typically desirable for the main parameter of interest. The nonparametic prior is 
found by an iterative Empirical Bayes type method, which aims to maximize marginal likelihood. The nonparametric prior 
may contain a point mass (default)}
\usage{
NonParaUpdatePrior(fitall,fitall0=NULL,  modus = "fixed", shrinkpara = NULL, shrinklc = NULL, lincombs=NULL, 
includeP0 = TRUE, unimodal = TRUE, logconcave = FALSE, symmetric = FALSE, allow2modes=TRUE, maxiter = 10, tol = 0.005, 
ndrawtot = 10000, ntotals = c(500, 2000, 5000, 10000), priornew = NULL, gaussinit0 = NULL, 
gammainit0 = NULL, p0init=0.8, p0lower = 0.5, pointmass = 0, maxsupport = 6, plotdens = TRUE, 
ncpus = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fitall}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} or from \code{\link[ShrinkBayes]{CombinePosteriors}}
}

\item{fitall0}{
An optional 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}} containing the fits under the null-model.
}

  \item{modus}{
Character string. Either \code{"fixed"}, \code{"random"} or \code{"logdisp"}. Type of variable for which the nonparametric prior is fit. 
}
  \item{shrinkpara}{
Character or character vector. Name(s) of the variable(s) for which the nonparametric prior is fit. Corresponding variable can be a factor. }

\item{shrinklc}{
Character string. Name of the linear combination for which the nonparametric prior is fit.
}

\item{lincombs}{
List object. Name of the list object that contains the linear combination(s), usually created by \code{inla.make.lincomb} or \code{AllComp}. 
Only required when \code{fitall} does not consist of INLA outputs, e.g. after \code{CombinePosteriors}.
}

 
\item{includeP0}{
Boolean. Should the nonparametric prior contain a point mass?
}
 
  \item{unimodal}{
Boolean. Should the nonparametric prior be unimodal?
}
  \item{logconcave}{
Boolean. Should the nonparametric prior be log-concave?
}
  \item{symmetric}{
Boolean. Should the nonparametric prior be symmetric?
}

\item{allow2modes}{
Boolean. Can the smooth part of the nonparametric prior contain a mode on the negative and positive half-plane?
}
  \item{maxiter}{
Integer. The maximum number of iterations in the iterative marginal procedure
} 
  \item{tol}{
Numeric. Tolerance for the Kolmogorov-Smirnov distance between two succesive iterations.
If this distance is smaller than \code{tol} the algorithm stops.
}
  \item{ndrawtot}{
Integer. Total number of draws from the empirical mixture of posteriors on which the new prior is fitted.
}
  \item{ntotals}{
Vector of integers. Consecutive number of posteriors that are used to determine the new prior. 
}
  \item{priornew}{
Numeric matrix or dataframe with 2 columns. Starting prior. Need not te be specified.
}
  \item{gaussinit0}{
Numeric vector of length 2. Mean and precision of Gaussian starting prior. Need not to be specified.
}
  \item{gammainit0}{
Numeric vector of length 2. Shape and rate of the Gamma starting prior. 
 Need not to be specified.
}

\item{p0init}{
Numeric between 0 and 1. Initial value of the mixture proportion of the point mass. 
Only relevant when \code{includeP0=TRUE}.
}
\item{p0lower}{
Numeric between 0 and 1. Lower bound for the mixture proportion of the point mass. Only relevant when \code{includeP0=TRUE}.
}  

\item{pointmass}{
Numeric. Location of the pointmass. Only relevant when \code{includeP0=TRUE}.
}

  \item{maxsupport}{
Numeric. maximum of the support of the prior and posteriors. For numerical stability. -maxsupport is the minimum of the support. 
}
  \item{plotdens}{
Boolean. Should the consecutive nonparametric estimates be plotted on the screen?
}
  \item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}
}

\details{
About \code{shrinklc}: it is assumed that only one type of linear combinations is present in the fit object \code{fitall}. 
Fitting a \code{unimodal} prior requires the package \code{"Iso"}, whereas a log-concave prior requires the package \code{"logcondens"} 
(Lutz and Rufibach, 2011). These options may be used to increase stability of the prior, in particular for the tails.}
\value{
A list with the following components:
 \item{priornew}{The (smooth) part of the nonparametric prior}
 \item{p0est}{Estimate of the point mass mixture component}
 \item{inputpar}{List of input parameters}
 \item{densall}{Consecutive fits of the nonparametric prior}
 \item{ksmaxall}{Kolmogorov-Smirnov distances between consecutive fits}
 \item{sdall}{Standard deviations for consecutive fits}
 \item{totalloglikall}{Marginal log-likelihood sums for consecutive fits}
 \item{nallforlik}{Number of elements used for computing \code{totalloglikall}}
 \item{notfit}{IDs of list elements of \code{fitall} for which a fit is not available}
}
\references{
Lutz, D. and Rufibach, K. (2011). logcondens: Computations related to univariate log-concave
density estimation. J. Statist. Software 39, 1-28.

Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}
\author{
Mark A. van de Wiel
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{MixtureUpdatePrior}} for parametric mixture priors (with point mass) and \code{\link{NonParaUpdatePosterior}} 
for computing posteriors from the output of this function. 
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}

