\name{AllComp}
\alias{AllComp}
\title{
Creating a multiple comparison object for use in INLA
}
\description{
This function creates all comparisons that do not include the baseline level from a factor object with 3 or more levels or from 
a character vector containing names of variables for which contrasts are desired.
}
\usage{
AllComp(vars)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vars}{A factor object with three or more levels or a character vector with names of model variables}
}
\details{
This is a convenience function to avoid having to specify each comparison separately. It calls 
\code{\link[INLA]{inla.make.lincomb}}. Please first use \code{\link[ShrinkBayes]{BaselineDef}} to define the baseline level. 
Comparisons with the baseline are default present in the output and hence do not need to be included.
}
\value{
A list object of the required format for the \code{lincomb} option in \code{\link[ShrinkBayes]{FitAllShrink}}, \code{\link[ShrinkBayes]{FitInlaAll}} and 
 \code{\link[INLA]{inla}}, 
}

\author{
Mark A. van de Wiel
}


\seealso{
\code{\link[INLA]{inla.make.lincomb}}, \code{\link[ShrinkBayes]{BaselineDef}}
}

