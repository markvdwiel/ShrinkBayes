\name{ShrinkBayesWrap}
\alias{ShrinkBayesWrap}
\title{
ShrinkBayes wrapper
}
\description{
This function enables applying ShrinkBayes using a one-line command. It includes the \code{ShrinkSeq()} (counts) or \code{ShrinkGauss()} (Gaussian data) functions,
the \code{FitAllShrink()} function, the \code{MixtureUpdatePrior()} function, the \code{MixtureUpdatePosterior()} function and the 
\code{SummaryTable()} function.
}
\usage{
ShrinkBayesWrap(data, form, paramtotest = NULL, allcontrasts = FALSE, multvscontrol = FALSE, fams = NULL, notfit = NULL, ncpus2use = 8, priorsimple = FALSE, approx0 = TRUE, diffthr = 0, direction = "two-sided", saveposteriors = TRUE, fileposteriors = "Posteriors.Rdata", ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{
Matrix or data frame containing the data. Rows are features, columns are samples. 
}
  \item{form}{
Formula starting with \code{~}. See \code{\link[INLA]{inla}} and 
\code{\link[INLA]{f}} for specification of the model formula.
}
  \item{paramtotest}{
Character. Name of variable to test. If \code{NULL} defaults to the first variable in \code{form}.
}
  \item{allcontrasts}{
Boolean. If \code{TRUE} all contrasts are tested if \code{paramtotest} is a factor with more than 2 levels.
}
  \item{multvscontrol}{
Boolean. If \code{TRUE} all comparisons with the first level of \code{paramtotest} are tested.
}
  \item{fams}{
Character string. Only relevant for count data. Either equal to \code{"poisson"}, \code{"zip"} 
(zero-inflated Poisson), \code{"nb"} (negative binomial), or the \code{"zinb"}
(zero-inflated negative binomial): likelihood to be used. Defaults to \code{"zinb"}.
}
  \item{notfit}{
Numeric vector of integers. Indices of data rows for which fitting is NOT required. 
}
  \item{ncpus2use}{
Integer. Number of cpus to use for calculations.
}
  \item{priorsimple}{
Boolean. If \code{TRUE} a simple pointmass-Normal prior is fit instead of the default a mixture. 
}
  \item{approx0}{
Boolean. If \code{TRUE} Savage-Dickey approximation is used for the null model (faster).
}
  \item{diffthr}{
Numeric. Threshold for the null-probability to compute. It demarcates the null-domain from the alternative domain.
}
  \item{direction}{
Character string. Use \code{direction="greater"} if the null-hypothesis if of the form \dQuote{parameter >= thr},
\code{direction="lesser"} if the null-hypothesis if of the form \dQuote{parameter <= thr}, use \code{direction="two-sided"}
if the null-hypothesis is of the form \dQuote{parameter = thr}, use \code{direction="equal"} for multi-parameter inference.
}
  \item{saveposteriors}{
Boolean. If \code{TRUE} posterior distributions of the parameters of interest are saved.
}
  \item{fileposteriors}{
Character. Name of file containing the posterior distributions
}
  \item{\dots}{
Further arguments passed on to \code{\link[ShrinkBayes]{ShrinkGauss}} (Gaussian data) or to \code{\link[ShrinkBayes]{ShrinkSeq}} (Counts).
See details. 
}
}
\details{
About \code{notfit}: note that this is not the same as removing the corresponding rows from the data. For correct fitting of the priors
it is essential to use all data. For final fitting, however, computing time can be considerably reduced if a quick screen was used, e.g. by 
\code{\link[ShrinkBayes]{ScreenData}}. For \code{saveposteriors}: please note that he object containing the posterior distributions 
can be large. Note that \code{direction} is automatically set to \code{equal} when multi-parameter inference is performed.
Specification of further arguments is particularly useful for shrinking other parameters of the model than \code{paramtotest} and the 
dispersion-related parameters (which are shrunken by default). In addition, the \code{ntag} and \code{maxiter} arguments (see 
\code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}) can be used to reduce computing time when testing.
}

\value{
List containing:
\item{FDRs}{Matrix of (B)FDRs, lfdrs, and effect size estimates}
\item{nsigsFDR01}{Numeric containing the number of significant results for each comparison at (B)FDR<=0.1}
\item{shrink}{Output of \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}} function}
\item{prior}{Output of \code{\link[ShrinkBayes]{MixtureUpdatePrior}} function, if applicable, \code{NULL} otherwise}
}

\references{
Van de Wiel MA, Neerincx M, Buffart TE, Sie D, Verheul HMW (2014).
ShrinkBayes: a versatile R-package for analysis of count-based sequencing data in complex study designs. 
BMC Bioinformatics. 15(1):116.

Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics. 14, 113-128.
}
\author{
Mark A. van de Wiel
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[ShrinkBayes]{ShrinkSeq}}, \code{\link[ShrinkBayes]{ShrinkGauss}}, \code{\link[ShrinkBayes]{MixtureUpdatePrior}},
\code{\link[ShrinkBayes]{ScreenData}}, \code{\link[ShrinkBayes]{BFUpdatePosterior}}, \code{\link[ShrinkBayes]{SummaryTable}},
\code{\link[ShrinkBayes]{plotPoster}}
}
\examples{
\dontrun{
#To empty current R memory and enlarge memory.limit (Windows) consider, before applying 'ShrinkBayes', using 
#rm(list=ls());gc();
#memory.limit(size = 4000)

#FIRST EXAMPLE
library(ShrinkBayes)
data(CAGEdata10000)
CAGEdata <- CAGEdata10000

#FOR FIRST TRY: USE 500 DATA ROWS ONLY
CAGEdata <- CAGEdata[1:500,]

#loads design of the brain study; includes 
data(design_brain)

#retrieves covariates from the design matrix
pers <- design_brain$pers  #persons
batch<-design_brain$batch   #batch
groupfac <- design_brain$groupfac #group (= brain region)

#Formula as used by INLA
form = ~ 1 + groupfac + batch + f(pers,model="iid") #'batch' and 'group' as fixed effect, 'pers' as random effect

#defaults to K-sample testing
SBcage <- ShrinkBayesWrap(CAGEdata,form)  

#Multiple Comparisons:
SBcagemc <- ShrinkBayesWrap(CAGEdata,form, allcontrasts = TRUE)  


#2ND EXAMPLE 
library(ShrinkBayes)
data(mirseqnorm)
data(designmirseq)

#Then, retrieve the covariates:
PM <- designmirseq$PM
indiv <- designmirseq$indiv
timepos <- designmirseq$timepos
chemo <- designmirseq$chemo
organ <- designmirseq$organ

#FORMULA; indiv IS A RANDOM EFFECT
form = ~ 1 + PM + timepos + chemo + organ + f(indiv)

#TRY IT FOR FIRST 200 DATA ROWS. USE ADDITIONAL SHRINKAGE;
#LIMIT ntag AND maxiter to save computing time; USE DEFAULTS FOR FINAL RUN
SBmir <- ShrinkBayesWrap(mirseqnorm[1:200,],form,shrinkaddfixed=c("organ","chemo","timepos"),ntag=c(50,100),maxiter=5,priorsimple=TRUE)  #2-sample testing for PM, with additional shrinkage for other fixed effects
}
}

