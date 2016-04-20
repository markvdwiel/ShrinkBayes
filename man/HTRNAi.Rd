\name{HTRNAi}
\alias{HTRNAi}
\docType{data}
\title{
A data set with 960 rows and 6 samples containing normalized High-Throughput RNAi data (Gaussian)}
\description{
This data set is used in the Van de Wiel et al. paper (see below) for illustrations.
}
\usage{data(HTRNAi)}
\format{
  The format is:
 num [1:960, 1:6] 0.54317601  0.19539894  0.38832444 -0.03428689  0.60492070 -0.17908612 
 - attr(*, "dimnames")=List of 2
  ..$ : chr [1:10000] "siRNA1" "siRNA2" "siRNA3" "siRNA4" ...
  ..$ : chr [1:6] "s1" "s2" "s3" "s4" "s5" "s6"
}
\details{
This data set is used to study sensitization of treatment effect by gene silencing in a cancer cell line using cell viability as read-out. 
The features represent (pools of) siRNA, where each pool targets the same mRNA. The data is assumed to be Gaussian distributed; 
hence, the examples below can also be used for typical mRNA/miRNA microarray studies. See the example in \code{\link{ShrinkGauss}} for the 
design of this study.
}

\references{
Van de Wiel MA,  De Menezes RX, Siebring-van Olst E, Van Beusechem VW 
(2013). Analysis of small-sample clinical genomics studies 
using multi-parameter shrinkage: application to high-throughput 
RNA interference screening. BMC Med Genom. 
}

\examples{
data(HTRNAi)
}
\keyword{datasets}
