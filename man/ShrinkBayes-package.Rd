\name{ShrinkBayes-package}
\alias{ShrinkBayes-package}
\alias{ShrinkBayes}
\docType{package}
\title{
Bayesian analysis of high-dimensional omics data, either Gaussian or counts 
}
\description{
This package provides differential expression analysis for a variety of omics data, including:
RNA-seq, CAGE, miRNA-seq, HT RNAi, mRNA microarray, miRNA microarray. It is particulary useful for relatively small 
sample sizes, because it applies Empirical Bayes-type estimation of priors for 
multi-parameter shrinkage to increase power and reproducibility. Mixture and nonparametric priors are accommodated.  
In addition, it provides 
Bayesian multiplicity correction, by lfdr and BFDR estimation. 
}
\details{
\tabular{ll}{
Package: \tab ShrinkBayes\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2012-10-15\cr
License: \tab GPL \cr
}

Functions: \cr
\code{\link{AllComp}}: Creates all comparisons that do not include the baseline level from a factor object with 3 or more levels \cr
\code{\link{BaselineDef}}: Defines the baseline level of factor variable. \cr
\code{\link{BFDR}}: BFDR computation, one-sided, two-sided and multiple comparisons \cr
\code{\link{BFUpdatePosterior}}:Updates posteriors for testing nested models that differ by more than one variable or simply extracts all posteriors for one 
particular parameter \cr
\code{\link{CombinePosteriors}}: Combines posteriors from two models for a given parameter and, possibly, for related linear combinations.\cr
\code{\link{CreateBlocks}}: Creates blocks of consecutives row indices from the row-dimension of a data set
and a block size \cr
\code{\link{FastScreenP}}: Screens and filters the data set using the data set \cr
\code{\link{FitAllShrink}}: Applies \code{\link[INLA]{inla}} to multiple data rows using shrinkage priors. \cr
\code{\link{MixtureUpdatePrior}}: Replacing a one-component prior by an optimal mixture prior \cr
\code{\link{MixtureUpdatePosterior}}: Updates posteriors from a given mixture prior \cr
\code{\link{NonParaUpdatePrior}}: Replacing a parametric prior by an optimal non-parametric one \cr
\code{\link{NonParaUpdatePosterior}}: Updates posteriors from a given non-parametric prior \cr
\code{\link{plotPoster}}: Plots posterior densities \cr
\code{\link{ReDefMiss}}: Redefines missing variables \cr
\code{\link{ScreenData}}: Fast p-value computation for simple 2-sample or K-sample designs. To be used as an initial screen. \cr
\code{\link{ShrinkBayesWrap}}: ShrinkBayes wrapper, apply ShrinkBayes using a one-line command \cr
\code{\link{ShrinkSeq}}: Simultaneous shrinkage by empirically estimating multiple priors,
Poisson and Negative Binomial likelihoods plus zero-inflated versions thereof
\cr
\code{\link{ShrinkGauss}}: Simultaneous shrinkage by empirically estimating multiple priors,
Gaussian likelihoods \cr
\code{\link{SummaryTable}}: Convenience function that produces a summary table \cr
\code{\link{SummaryWrap}}: Convenience function to compute either posterior null-probabilities (lfdr) or posterior means \cr
\cr
Data:\cr
\code{\link{CAGEdata10000}}: A data set with 10,000 rows containing normalized counts \cr
\code{\link{datsim}}: A simulated data set (Gaussian) with 1500 rows and 8 columns \cr
\code{\link{HTRNAi}}: A data set with 960 rows and 6 samples containing normalized High-Throughput RNAi data (Gaussian) \cr
\code{\link{mirseqnorm}}: A data set with 2060 rows and 55 samples containing normalized counts \cr
}

\author{
Mark A. van de Wiel, mark.vdwiel@vumc.nl

}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.

Van de Wiel MA,  De Menezes RX, Siebring-van Olst E, Van Beusechem VW 
(2013). Analysis of small-sample clinical genomics studies 
using multi-parameter shrinkage: application to high-throughput 
RNA interference screening. BMC Med Genom. 

Van de Wiel MA, Neerincx M, Buffart TE, Sie D, Verheul HMW (2014).
ShrinkBayes: a versatile R-package for analysis of count-based sequencing data in complex study designs. 
BMC Bioinformatics. 15(1):116.
}

