\name{CreateBlocks}
\alias{CreateBlocks}
\title{A function that creates blocks of consecutives row indices from the row-dimension of a data set
and a block size}
\description{This is a convenience function, mainly used for sequential running on very large data sets. 

}
\usage{
CreateBlocks(ds, size = 10000)
}
\arguments{
  \item{ds}{Matrix or data frame.          
}
  \item{size}{Numeric indicating the (maximum) size of the blocks
}
}
\details{
The last block may consist of less indices than \code{size}.
}
\value{A matrix, the rows of which present the starting and end indices of the blocks.
}

\author{
Mark A. van de Wiel
}

\examples{
\dontrun{See the example in the vignette
on large data sets.}
}

