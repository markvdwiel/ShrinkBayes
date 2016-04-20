\name{BaselineDef}
\alias{BaselineDef}
\title{Defines the baseline level of factor variable. 
}
\description{
Baseline level will be recoded to \code{0} to ensure 
that all further analyses correctly identify this level to be the baseline.
}
\usage{
BaselineDef(fact, baselinegroup = as.character(levels(get(fact)))[1])
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fact}{Factor variable}
  \item{baselinegroup}{Character string defining the baseline level. Should be a level of \code{fact}.
}
}
\details{
If a level \code{0} already exists, this will be redefined into level \code{0b}.
}
\value{A factor variable with level 0 as the baseline
}

\author{
Mark A. van de Wiel
}


\seealso{
\code{\link[ShrinkBayes]{AllComp}}
}

