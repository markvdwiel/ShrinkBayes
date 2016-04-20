\name{CAGEdata10000}
\alias{CAGEdata10000}
\docType{data}
\title{
A data set with 10,000 rows containing normalized counts}
\description{
This data set is used in the Van de Wiel et al. paper (see below) for illustrations.
}
\usage{data(CAGEdata10000)}
\format{
  The format is:
 num [1:10000, 1:25] 0 0 0 0 2 8 35 2 31 0 ...
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:10000] "Tag.1" "Tag.2" "Tag.3" "Tag.4" ...
  ..$ : chr [1:25] "raw.0325_Frontal" "raw.0325_Hippocampus" "raw.0325_Putamen" "raw.0325_Temporal" ...
}
\details{
This data set is part of a larger data set. It is used to study promotor activity. The features represent tag clusters. 
It corresponds to the design \code{\link{design_brain}}.
}
\source{
Pardo LM, Rizzu P, Francescatto M, Vitezic M, Leday GG, Sanchez JS, Khamis A, Takahashi H, van de Berg WD, Medvedeva YA, van de Wiel MA, Daub CO, Carninci P, Heutink P (2013). 
Regional differences in gene expression and promoter usage in aged human brains. Neurobiol Aging. 34(7):1825-36
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}
\examples{
data(CAGEdata10000)
data(design_brain)
## maybe str(CAGEdata10000) ; plot(CAGEdata10000) ...
}
\keyword{datasets}
