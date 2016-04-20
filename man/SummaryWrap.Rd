\name{SummaryWrap}
\alias{SummaryWrap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience function to compute either posterior null-probabilities (lfdr) or posterior means}
\description{
This function allows one to summarize the final posteriors. It also serves as input for the computation of Bayesian False Discovery Rate (
see \code{\link{BFDR}}). 
}
\usage{
SummaryWrap(posteriorlist, thr = 0, ncpus = 1, updateby = 10000, summary = "lfdr", direction = "two-sided", pointmass=0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{posteriorlist}{
A list of lists where each sublist contains posteriors of (possibly) multiple parameters (or linear combinations) from the same data. 
A posterior is represented as a matrix with two columns denoting support and posterior density.
}
  \item{thr}{
Numeric. Threshold \code{thr} for the null-probability to compute. It demarcates the null-domain from the alternative domain. Only relevant when 
\code{summary="lfdr"}
}
  \item{ncpus}{Integer.  The number of cpus to use for parallel computations.
}
  \item{updateby}{Integer. Show an update for each \code{updateby} number of posteriors executed.
}
  \item{summary}{
Character string. Use \code{summary="lfdr"} to obtain posterior null-probabilities (lfdrs), \code{summary="postmean"} to compute
posterior means and \code{summary="postpi0"} to obtain (joint) posterior point null-probability 
(may equal lfdr when \code{direction="equal"} is used)
}
\item{direction}{
Character string. Use \code{direction="greater"} if the null-hypothesis if of the form \dQuote{parameter >= thr},
\code{direction="lesser"} if the null-hypothesis if of the form \dQuote{parameter <= thr}, use \code{direction="two-sided"}
if the null-hypothesis is of the form \dQuote{parameter = thr}, use \code{direction="equal"} for multi-parameter inference.
}
\item{pointmass}{Numeric. Location of the pointmass. Only relevant when the components of posteriorlist contain a point mass}
}
\details{
About \code{posteriorlist}: this usually results from \code{\link{NonParaUpdatePosterior}}, \code{\link{MixtureUpdatePosterior}}, 
\code{\link{BFUpdatePosterior}} or \code{\link{SpikeSlabUpdatePosterior}}. 
It may either contain only continous posteriors or a mix of those with a point mass (usually at 0 to force selection properties). 
Please note that posterior null-probabilities can be interpreted as (a Bayesian version) of local false discovery rates (lfdr), hence these are obtained 
when setting \code{summary="lfdr"}. These lfdrs serve as inpute for computing Bayesian False Discovery Rates, see \code{\link{BFDR}}. 

We recommend to use \code{direction="equal"} only for Bayes factor type inference 
(for comparing nested models that differ by more than one parameter). If point null inference is desired for a single parameter, we 
recommend to use both \code{direction="greater"} and \code{direction="lesser"} as outlined in the examples.

}
\value{
Numerical matrix. The number of rows equals the length of \code{posteriorlist}, the number of columns equals the number of posteriors in each 
sublist.  
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.

Efron B, Tibshirani R, Storey JD and Tusher V. (2001). Empirical bayes analysis
of a microarray experiment. J. Amer. Statist. Assoc. 96, 1151-1160.

}
\author{
Mark A. van de Wiel}
\note{
Functionality may be extended in the future}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SummaryTable}}, \code{\link{NonParaUpdatePosterior}}, \code{\link{MixtureUpdatePosterior}} and \code{\link{BFDR}}.
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}

