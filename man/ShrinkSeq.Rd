\name{ShrinkSeq}
\alias{ShrinkSeq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simultaneous shrinkage by empirically estimating multiple priors,
Poisson and Negative Binomial likelihoods plus zero-inflated versions thereof
}
\description{
This function implements the iterative joint procedure as described in Van
de Wiel et al. (2012). It recursively applies INLA, and fits parameters
of multiple priors by an EM-type algorithm.
}
\usage{
ShrinkSeq(form, dat, fams = "zinb", shrinkfixed = NULL, 
shrinkaddfixed = NULL, shrinkrandom = NULL, shrinkaddrandom=NULL,shrinkdisp = TRUE, 
shrinkp0 = FALSE, mixtdisp = FALSE, mixtrand = FALSE, excludefornull=NULL, fixedmeanzero = FALSE, 
addfixedmeanzero = TRUE, curvedisp = FALSE, maxiter = 15,  
ntag=ifelse(is.null(excludefornull),c(100,200,500,1000),c(1000)), ntagcurve = 5000, 
fixed = c(0, 1/10), addfixed = c(0, 1/10), randomprec = c(1,10^(-5)), addrandomprec = c(1,10^(-5)),
logdisp = c(0, 0.01), diracprob0=ifelse((mixtrand | !is.null(excludefornull)), 0.8, 0.2), logitp0 = c(0, 0.01), 
fixedseed = TRUE, ndraw = 10000, safemode = TRUE, tol = ifelse((mixtrand | !is.null(excludefornull)),0.005,0.01), 
tolrand = 0.02, mliktol = 0.1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{form}{
Formula. See \code{\link[INLA]{inla}} and 
\code{\link[INLA]{f}} for specification of the model formula.
}
  \item{dat}{
Matrix or data frame containing the data. Rows are features, columns are samples. 

}
  \item{fams}{
Character string. Either equal to \code{"poisson"}, \code{"zip"} 
(zero-inflated Poisson), \code{"nb"} (negative binomial), or \code{"zinb"}
(zero-inflated negative binomial): likelihood to be used.
}
  \item{shrinkfixed}{
String. Name of the fixed parameter for which a Gaussian prior 
is fit. } 

  \item{shrinkaddfixed}{
String or character vector. Name(s) of additional fixed parameter(s) for which a Gaussian  
prior is fit. } 

  \item{shrinkrandom}{
String. Name of Gaussian random effects parameter for which an 
inverse-Gamma prior is fit to the variance} 

\item{shrinkaddrandom}{
String or character vector. Name(s) of additional random effects parameter(s) for which 
inverse-Gamma priors are fit to the variances} 

  \item{shrinkdisp}{
Boolean. Fit a prior to the log-overdispersion parameter of the 
(zero-inflated) negative binomial?} 

  \item{shrinkp0}{
Boolean. Fit a Gaussian prior to the logit(zero inflation) 
parameter of the zero-inflated poisson or negative binomial? } 

  \item{mixtdisp}{
Boolean. Add a point mass on zero to the Gaussian prior
  of the log-overdispersion parameter? Should be FALSE when \code{mixtrand=TRUE}.
}
  \item{mixtrand}{Boolean. Add a point mass on zero to the inverse-Gamma prior
  for the random effect? Should be FALSE when \code{mixtdisp=TRUE}.
}

\item{excludefornull}{
String or character vector. Name(s) of parameter(s) to be excluded from the given model in \code{form} resulting in the null mode
}

  \item{fixedmeanzero}{
Boolean. Fix the mean of the Gaussian prior of the 
\code{shrinkfixed} parameter to zero?} 

  \item{addfixedmeanzero}{
Boolean or list of booleans. Fix the mean(s) of the Gaussian prior of the 
\code{shrinkaddfixed} parameter(s) to zero?} 

  \item{curvedisp}{
Boolean. Apply a curvature prior for log-overdispersion of the 
negative binomial? 
 } 
 
  \item{maxiter}{
Integer. The maximum number of iterations in the iterative marginal procedure}
  
  \item{ntag}{
Integer vector. Consecutive number of data rows used for fitting the priors. Rows are sampled at random with a fixed seed. 
}
  \item{ntagcurve}{
Integer vector. Only relevant if \code{curvedisp=TRUE}. Consecutive number of data rows used for fitting a curvature prior for overdispersion. 
Rows are sampled at random with a fixed seed. 
}
   \item{fixed}{
Numeric vector of length 2. Mean and precision of initial 
Gaussian prior of \code{shrinkfixed}} 

  \item{addfixed}{
Numeric vector of length 2 or list of vectors of size 2. Mean(s) and precision(s) of initial 
Gaussian prior(s) of \code{shrinkaddfixed}}

  \item{randomprec}{
Numeric vector of length 2. Shape and rate of initial 
inverse-Gamma prior of \code{shrinkrandom}} 

\item{addrandomprec}{
Numeric vector of length 2. Shape and rate of initial 
inverse-Gamma priors of \code{shrinkaddrandom}}

  \item{logdisp}{
Numeric vector of length 2. Mean and precision of initial 
Gaussian prior of log-overdispersion} 
 
  \item{diracprob0}{
Numeric between 0 and 1. Initial mixture proportion for 
zero-component of the prior for \code{shrinkdisp} or for the 
random effects variance. Only relevant when 
\code{mixtdisp=TRUE} or \code{mixtrand=TRUE}} 

  \item{logitp0}{
Numeric vector of length 2. Mean and precision of initial 
Gaussian prior of logit (zero-inflation) parameter}
 
  \item{fixedseed}{
Boolean. Use a fixed seed for sampling the features?}

  \item{ndraw}{
Integer. The number of random draws from the empirical mixture 
of posteriors } 

  \item{safemode}{
Boolean. If \code{TRUE}, only use features for which both models 
provide a fit, in case of a mixture prior. 
} 

\item{tol}{
Numeric. Tolerance for the Kolmogorov-Smirnov distance between two succesive iterations.
}
  \item{tolrand}{
Numeric. Tolerance for the Kolmogorov-Smirnov distance between two succesive iterations for the prior of the random effect variance .
}

\item{mliktol}{
Numeric. Relative tolerance in percentage between two succesive iterations for the total marginal likelihood. 
}

  \item{\dots}{
Further arguments passed on to \code{\link[INLA]{inla}}. 
}
} 


\details{If a prior is fit for a particular parameter, this 
effectuates Bayesian shrinkage for that parameter. The 
iterative algorithm operates on an increasing number of tags, 
which are indicated in \code{ntag}. The algorithm stops if the 
following conditions are satisfied: 1. it has processed on the 
largest element of \code{ntag} 2a. It has reached 
\code{maxiter} OR 2b. For all non-random effects priors the 
Kolmogorov-Smirnov (KS) distance between two consecutive fits 
of the prior is smaller than \code{tol} AND the KS distance 
between two consecutive fits of the random effects variance is 
smaller than \code{tolrand} OR 2c. The relative difference (Total log marg lik - Previous total log marg lik)/
abs(Total log marg lik) < mliktol. 

Computing time is largely 
determined by the largest element in \code{ntag} and 
\code{maxiter}. We advise to use at least \code{max(ntag)=1000} 
and \code{max(ntag)=2500} when either \code{mistdisp=TRUE} or 
\code{mixtrand=TRUE}. About \code{shrinkfixed} and 
\code{shrinkaddfixed}: if one of the fixed parameters is the 
primary parameter of interest use \code{shrinkfixed} for this 
one, and, when desired, \code{shrinkaddfixed} for an another 
one.}


\value{
A list object containing:
  \item{pmlist}{Final estimates of the parameters of the (fitted) priors}
  \item{ksall}{Kolmogorv-Smirnov distances between consecutive estimated priors}
  \item{paraall}{List object containing all consecutive parameter estimates}
  \item{curvedispfun}{Function that links log-overdispersion with log-mean; only relevant when \code{curvedisp=TRUE}}
  \item{inputpar}{Input parameters used}
  \item{addfixed}{Parameter values of additional fixed parameters. Passed on for convenience.}
  \item{addrandom}{Parameter values of additional random effects parameters. Passed on for convenience.}
  \item{typelik}{Type of likelihood used}
}

\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.

Van de Wiel MA, Neerincx M, Buffart TE, Sie D, Verheul HMW (2014).
ShrinkBayes: a versatile R-package for analysis of count-based sequencing data in complex study designs. 
BMC Bioinformatics. 15(1):116.
}

\author{
Mark A. van de Wiel}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link[ShrinkBayes]{ShrinkGauss}}
}

\examples{
\dontrun{
#For an example that combines posteriors (using a mixture prior for overdispersion like in the CAGE example in Van de Wiel et al., 2012): 
#see CombinePosteriors
#To empty current R memory and enlarge memory.limit (Windows) consider, before applying 'ShrinkBayes', using 
#rm(list=ls());gc();
#memory.limit(size = 4000)

#loads brain sequencing data (CAGE); 10,000 rows (tag clusters), 25 samples
data(CAGEdata10000)
CAGEdata <- CAGEdata10000

#FOR FIRST TRY: USE 1000 DATA ROWS ONLY
CAGEdata <- CAGEdata[1:1000,]

#loads design of the brain study; includes covariates
data(design_brain)

#retrieves covariates from the design matrix
pers <- design_brain$pers  #persons
batch<-design_brain$batch   #batch
groupfac <- design_brain$groupfac #group (= brain region)

#Number of cpus to use in (parallel) computations
ncpus2use <- 10

#redefines baseline level as '0' such that it is in further analysis always used as baseline. 
groupfac <- BaselineDef("groupfac",baselinegroup="1")

#By default only comparisons with the baseline are included. Use the following convenience function to create the other 
#pair-wise comparisons for a given factor when this contains 3 or more levels (groups)
lincombvec <- AllComp("groupfac")

#Formula as used by INLA
form = y ~ 1 + groupfac + batch + f(pers,model="iid") #'batch' and 'group' as fixed effect, 'pers' as random effect

#Simultaneous hrinkage for 'groupfac' and 'pers'. In addition, the negative binomial overdispersion is shrunken by default (see shrinkdisp argument). 
#In the example below we allow a mixture prior for overdispersion (mixtdisp=TRUE). This increases computing time by a factor 2.
#INCREASE maxiter FOR OBTAINING FINAL RESULTS (DEFAULT = 15)
shrinksimul <- ShrinkSeq(form=form, dat=CAGEdata,shrinkfixed="groupfac",shrinkrandom="pers",ncpus=ncpus2use, maxiter=3)

#Fit all using the priors resulting from ShrinkSeq and zi-nb likelihood; also find posteriors for the linear combinations defined above
fitzinb <- FitAllShrink(form,dat=CAGEdata,fams="zinb",shrinksimul,ncpus=ncpus2use,lincomb=lincombvec)

#Find nonparametric prior for all group-wise differences. These include the linear combinations. 
#INCREASE maxiter FOR OBTAINING FINAL RESULTS (DEFAULT = 15)
npprior <- NonParaUpdatePrior(fitall=fitzinb,modus="fixed", includeP0 = FALSE, shrinkpara="groupfac", shrinklc=TRUE,ncpus=ncpus2use, maxiter = 3)

#Update posteriors using the nonparametric prior
nppostshr <- NonParaUpdatePosterior(fitzinb,npprior,ncpus=ncpus2use)

#Compute local fdrs ( = posterior null (tail-)probabilities under shrinkage prior)
#Here P(theta <= 0.1) is computed, assuming {theta<=0.1} is the null-domain.
lfdr <- SummaryWrap(nppostshr, thr = 0.1)

#Compute Bayesian FDRs; all BFDRs for each comparison [matrix format]
BFDRs <- BFDR(lfdr)

#Compute Bayesian FDRs; BFDR for multiple comparisons
BFDRmult <- BFDR(lfdr,multcomp=TRUE)
}
}




