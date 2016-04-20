\name{datsim}
\alias{datsim}
\docType{data}
\title{
A simulated data set with 1500 rows and 8 columns}
\description{
This simulated data set is used to illustrate several functions in this package
}
\usage{data(CAGEdata10000)}
\format{
  The format is:
 num [1:1500, 1:8] -0.04191917  0.070037754  0.8198383 -0.35443640 ...
}
\details{
For reproducibility reasons we make this simulated data set available. 
It was generated as follows: \cr
\code{Specifications: 1500 rows, 2*4 samples, fraction non-differential = 2/3 (first 1000 rows)} \cr 
\code{effect size distribution for differential features: N(0,1), noise sd: 0.5.} \cr
\code{datsim1 <- matrix(rnorm(8000,mean=0,sd=0.5),nrow=1000)} \cr
\code{meanvec <- matrix(rep(rnorm(500,0,1),4),nrow=500)} \cr
\code{datsim2 <- cbind(matrix(rnorm(2000,mean=0,sd=0.5),nrow=500),matrix(rnorm(2000,mean=0,sd=0.5),nrow=500) + meanvec)} \cr
\code{datsim <- rbind(datsim1,datsim2)} 
}

\seealso{\code{\link{HTRNAi}}, \code{\link{CAGEdata10000}}}

\examples{
data(datsim)
}
\keyword{datasets}
