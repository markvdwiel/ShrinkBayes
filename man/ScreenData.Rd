\name{ScreenData}
\alias{ScreenData}
\title{
Fast p-value computation for simple 2-sample or K-sample designs. 
}
\description{
This function provides fast p-value combination for very large data sets based on t-tests, F-tests, Wilcoxon two-sample tests or Kruskal-Wallis tests.  
It may be used as an initial screen and in combination with 
\code{\link{FastScreenP}} to filter data and speed up computations.
}
\usage{
ScreenData(dat, variable, np = TRUE, ncpus = 2, updateby = 10000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{dat}{
Matrix or data frame containing the data. Rows are features, columns are samples. 
}
  \item{variable}{
String. Name of the factor variable that contains the group information for each sample
}
  \item{np}{
Boolean. Should a nonparametric test be used or not?
}
  \item{ncpus}{
Integer.  The number of cpus to use for parallel computations.
}
  \item{updateby}{
Integer. Show an update for each \code{updateby} number of features executed.
}
}
\details{
If \code{np=TRUE} then a Wilcoxon two-sample test is performed when the number of levels of \code{variable} equals 2, and 
a Kruskal-Wallis test is performed when the number of levels exceeds 2. If \code{np=FALSE}, a two-sample t-test or F-test is used 
instead.
}
\value{
Numerical vector containing the p-values.
}

\author{
Mark A. van de Wiel
}

\seealso{
\code{\link{wilcox.test}}, \code{\link{kruskal.test}}, \code{\link{t.test}}, \code{\link{anova}}.
}

\examples{
# See FastScreenP
}

