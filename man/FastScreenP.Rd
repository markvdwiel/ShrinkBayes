\name{FastScreenP}
\alias{FastScreenP}
\title{
Using p-values to screen and filter the data set
}
\description{
This is a convenience function to speed up your computations when the data set is very large. 
 
}
\usage{
FastScreenP(pvalues, method = "padjust", threshold = ifelse(method == "pvalue", 0.1, 0.5), adjmethod = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{pvalues}{
Numerical vector of p-values.
}
  \item{method}{
Character string either equal to "padjust" or "pvalue"
}
  \item{threshold}{
Numeric. Upper threshold for reporting indices
}
  \item{adjmethod}{Character string.
P-value adjustment method in \code{\link{p.adjust}} 
}
}
\details{
This function is particularly convenient when you have access to a method that provides very fast computation of p-values. If such p-values are large,
you may opt to not compute posteriors for the corresponding features. Please be aware that such filtering should take place after all priors have been 
fit and fixed to avoid biases. See the vignette for an example. 
\code{\link[ShrinkBayes]{ScreenData}} provides computation of p-values for two-group and multi-group settings.
}
\value{Integer vector containing all indices with (adjusted) p-values below the threshold
}

\author{
Mark A. van de Wiel
}

\note{
Use with CARE. Do not use it to filter data for \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}, because this would 
introduce a bias.}


\seealso{
\code{\link[ShrinkBayes]{ScreenData}} 
}

\examples{
# Simulation adopted from limma. Simulate gaussian data for 1000 probes and 6 samples
# Samples are in two groups
# First fifty probes are differentially expressed in second group
# Std deviations vary with prior df=4
sd <- 0.3*sqrt(4/rchisq(1000,df=4))
y <- matrix(rnorm(1000*6,sd=sd),1000,6)
rownames(y) <- paste("Feature",1:1000)
y[1:50,4:6] <- y[1:50,4:6] + 2

group <- factor(c(1,1,1,2,2,2))

#performs t-test for all 1000 rows
pvals <- ScreenData(y, "group", np = FALSE, ncpus = 2)

#screens on the basis of FDR, FDR <= 0.5. Since this is an initial screen, it is best to be liberal.
whichin <- FastScreenP(pvals)

#run ShrinkSeq, ShrinkGauss and UpdatePrior functions on the ENTIRE data set:
#for these functions computing time does not strongly depend on the number of features
#For function FitAllShrink() use the selected set of features:
\dontrun{
fit <- FitAllShrink(form=form,dat=y[whichin,],...)
}
}



