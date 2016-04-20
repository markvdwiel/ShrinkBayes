\name{SummaryTable}
\alias{SummaryTable}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Convenience function that produces a summary table}
\description{
This function produces various summary statistics from the final posteriors. Only features with BFDR smaller than \code{BFDRthr} are 
selected.
}
\usage{
SummaryTable(poster,BFDRthr=0.1,diffthr = 0,direction="two-sided",pointmass=0,ndigit=3,ncpus=1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{poster}{
A list of lists where each sublist contains posteriors of (possibly) multiple parameters (or linear combinations) from the same data.
A posterior is represented as a matrix with two columns denoting support and posterior density.
}
  \item{BFDRthr}{
Numeric. Selection threshold for the BFDR. 
}

\item{diffthr}{
Numeric. Threshold for the null-probability to compute. It demarcates the null-domain from the alternative domain.
}
  
\item{direction}{
Character string. Use \code{direction="greater"} if the null-hypothesis if of the form \dQuote{parameter >= thr},
\code{direction="lesser"} if the null-hypothesis if of the form \dQuote{parameter <= thr}, use \code{direction="two-sided"}
if the null-hypothesis is of the form \dQuote{parameter = thr}, use \code{direction="equal"} for multi-parameter inference.
}

\item{pointmass}{Numeric. Location of the pointmass. Only relevant when the components of posteriorlist contain a point mass.}

\item{ndigit}{Integer.  The number of decimal digits to be used for the output.
}

\item{ncpus}{Integer.  The number of cpus to use for parallel computations.
}
}
\details{
About \code{posteriorlist}: this usually results from \code{\link{NonParaUpdatePosterior}}, \code{\link{MixtureUpdatePosterior}}, 
\code{\link{BFUpdatePosterior}} or \code{\link{SpikeSlabUpdatePosterior}}.
It may either contain only continous posteriors or a mix of those with a point mass (usually at 0 to force selection properties).

The summary table contains: poster null-probabilities, local fdr (equal posterior null-probabilities of \code{direction ="equal"} and 
\code{diffthr = 0}), Bayesian FDR and posterior means for all parameters or contrasts for which the posteriors are provided in 
\code{posteriorlist}.

We recommend to use \code{direction="equal"} only for Bayes factor type inference
(for comparing nested models that differ by more than one parameter). If point null inference is desired for a single parameter, we
recommend to use both \code{direction="greater"} and \code{direction="lesser"} as outlined in the examples.
}
\value{
Data frame. The number of rows equals the number of features with \dQuote{BFDR <= BFDRthr}. The first column contains the indices of 
the selected features, the other columns contain the summary statistics
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012).
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.

}
\author{
Mark A. van de Wiel}
%\note{
%Functionality may be extended in the future}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{SummaryWrap}}, \code{\link{BFDR}} and \code{\link{plotPoster}} for plotting posterior densities.
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}

