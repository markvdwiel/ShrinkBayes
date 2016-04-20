\name{designmirseq}
\alias{designmirseq}
\docType{data}
\title{
Design corresponding to the mirseqnorm data 
}
\description{
Contains the design of the \code{\link[ShrinkBayes]{mirseqnorm}} data
}
\usage{data(designmirseq)}
\format{
  A data frame with 55 observations on the following 8 variables.
  \describe{
    \item{\code{PM}}{a factor with levels \code{P} \code{M} }
    \item{\code{indiv}}{a factor with 27 levels}
    \item{\code{timepos}}{Binary}
    \item{\code{chemo}}{Binary}
    \item{\code{organ1}}{Binary}
    \item{\code{organ2}}{Binary}
    \item{\code{organ3}}{Binary}
    \item{\code{organ}}{Factor variable}
  }
}
\details{
For details on the design, please see the references below. Note that 
\code{organ} contains the same info as \code{organ1,organ2,organ3}, but then as factor variable
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
data(designmirseq)
## maybe str(CAGEdata10000) ; plot(CAGEdata10000) ...
}
\keyword{datasets}

