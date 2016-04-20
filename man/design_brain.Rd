\name{design_brain}
\alias{design_brain}
\docType{data}
\title{
Design corresponding to the CAGE10000 data 
}
\description{
Contains the design of the \code{\link[ShrinkBayes]{CAGEdata10000}} data
}
\usage{data(design_brain)}
\format{
  A data frame with 25 observations on the following 3 variables.
  \describe{
    \item{\code{pers}}{a numeric vector}
    \item{\code{batch}}{a factor with levels \code{0} \code{1}}
    \item{\code{groupfac}}{a factor with levels \code{1} \code{2} \code{3} \code{4} \code{5}}
  }
}
\details{
For details on the design, please see the references below.
}
\source{
Pardo et al. (2012). 
}
\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}
\examples{
data(CAGEdata10000)
data(design_brain)
## maybe str(design_brain) ; plot(design_brain) ...
}
\keyword{datasets}
