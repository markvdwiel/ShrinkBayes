\name{mirseqnorm}
\alias{mirseqnorm}
\docType{data}
\title{
A miRseq data set with 2,060 rows containing normalized counts}
\description{
This data set is used in the Van de Wiel et al. paper (see below) for illustrations.
}
\usage{data(mirseqnorm)}
\format{
Numerical data frame of size 2060 x 55.
}
\details{
This data set is part of a larger data set. It is used to study promotor activity. The features represent tag clusters. 
It corresponds to the design \code{\link{design_brain}}.
}
\source{
van de Wiel et al. (2013). Submitted.
}
\references{
Van de Wiel MA, Neerincx M, Buffart TE, Verheul HMW (2013). 
ShrinkBayes: a versatile R-package for analysis
of count-based sequencing data in complex study designs. Submitted.
}
\examples{
data(mirseqnorm)
data(designmirseq)
## maybe str(CAGEdata10000) ; plot(CAGEdata10000) ...
}
\keyword{datasets}
