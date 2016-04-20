\name{BFDR}
\alias{BFDR}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
BFDR computation, one-sided, two-sided and multiple comparisons
}
\description{
This function computes Bayesian FDRs from posterior tail probabilities or local fdrs P(H0|Y).
}
\usage{
BFDR(lfdr1,lfdr2=NULL,multcomp=FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{lfdr1}{Vector or matrix containing tail probabilities for each comparison or parameter}
  \item{lfdr2}{Vector or matrix containing tail probabilities (to the the other side than those in \code{lfdr1}) for each comparison or parameter}
  \item{multcomp}{Boolean, TRUE or FALSE indicating whether a multiple comparison is desired}
}
\details{
See Ventrucci et al. reference below for interpretation of BFDR, also in comparison with ordinary FDR. Posterior tail probabilities can be 
computed with the function \code{\link[ShrinkBayes]{SummaryWrap}}
}
\value{
A numerical vector or matrix. A matrix is returned when the input \code{lfdrs} are matrices and \code{multcomp=FALSE}.
}
\references{
Ventrucci, M, Scott, EM and Cocchi, D (2011). Multiple testing on standardized
mortality ratios: a Bayesian hierarchical model for FDR estimation. Biostatistics 12, 51-67.

Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}
\author{
Mark A. van de Wiel
}

\seealso{
\code{\link[ShrinkBayes]{SummaryWrap}}
}

\examples{
#See ShrinkSeq, ShrinkGauss and CombinePosteriors
}
