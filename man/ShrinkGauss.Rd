\name{ShrinkGauss}
\alias{ShrinkGauss}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Simultaneous shrinkage by empirically estimating multiple priors,
Gaussian likelihoods
}

\description{
This function implements the iterative joint procedure as described in Van
de Wiel et al. (2012) and adapted to Gaussian likelihoods in 
Van de Wiel et al. (2013). It recursively applies INLA, and fits parameters
of multiple priors by an EM-type algorithm.
}

\usage{
ShrinkGauss(form, dat, maxiter = 10, shrinkfixed = NULL, 
shrinkaddfixed = NULL, shrinkrandom = NULL, shrinkaddrandom=NULL, shrinksigma = TRUE, 
mixtrand = FALSE, excludefornull=NULL, fixedmeanzero = FALSE, addfixedmeanzero = TRUE, 
ntag=ifelse(is.null(excludefornull),c(100,200,500,1000),c(1000)), fixed = c(0, 1/10), addfixed = c(0, 1/10), 
randomprec = c(1,10^(-5)),addrandomprec = c(1,10^(-5)), precerr = c(1,10^(-5)), diracprob0=ifelse((mixtrand | !is.null(excludefornull)), 0.8, 0.2), fixedseed = TRUE, ndraw = 10000, 
safemode = TRUE, tol=ifelse((mixtrand | !is.null(excludefornull)),0.005,0.01), tolrand = 0.02,  mliktol = 0.1, designlist=NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{form}{
Formula. See \code{\link[INLA]{inla}} and 
\code{\link[INLA]{f}} for specification of the model formula.
}
  \item{dat}{
Matrix, data frame or list containing the data. Rows are features, columns are samples. For lists: 
each component represents a feature.
}
  \item{maxiter}{
Integer. The maximum number of iterations in the iterative marginal procedure}
  
  \item{ntag}{
Integer vector. Consecutive number of data rows used for fitting the priors. Rows are samples at random with a fixed seed. 
}
  \item{shrinkfixed}{
String. Name of the fixed parameter for which a Gaussian prior 
is fit. } 

  \item{shrinkaddfixed}{
String or character vector. Name(s) of additional fixed parameter(s) for which Gaussians  
priors are fit. } 
 
  \item{shrinkrandom}{
String. Name of Gaussian random effects parameter for which an 
inverse-Gamma prior is fit to the variance} 

 \item{shrinkaddrandom}{
String or character vector. Name(s) of additional random effects parameter(s) for which 
inverse-Gamma priors are fit to the variances} 

  \item{shrinksigma}{
Boolean. Fit a prior to the error variance? } 
 
  \item{mixtrand}{Boolean. Add a point mass on zero to the inverse-Gamma prior
  for the random effect?
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

  \item{precerr}{
Numeric vector of length 2. Shape and rate of initial 
inverse-Gamma prior of the error variance} 

  \item{diracprob0}{
Numeric between 0 and 1. Initial mixture proportion for 
zero-component of the prior for \code{shrinkrandom}. Only 
relevant when \code{mixtrand=TRUE}} 

  \item{fixedseed}{
Boolean. Use a fixed seed for sampling the features?}

  \item{ndraw}{
Integer. The number of random draws from the empirical mixture 
of posteriors } 

  \item{safemode}{
Boolean. If TRUE, only use features for which both models 
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

\item{designlist}{
List. Components are data frames containting the variables in \code{form}. Length of list should equal nr of features in \code{data}. If \code{NULL} design is assumed to be the same for all features (and inferred from \code{form})
}
  
  \item{\dots}{
Further arguments passed on to \code{\link[INLA]{inla}}. 
} 
}

\details{If a prior is fit for a particular parameter, this 
effectuates Bayesian shrinkage for that parameter. The 
iterative algorithm operates on an increasing number of 
features (e.g. tags), which are indicated in \code{ntag}. The 
algorithm stops if the following conditions are satisfied: 1. 
it has processed on the largest element of \code{ntag} 2a. It 
has reached \code{maxiter} or 2b. For all non-random effects 
priors the Kolmogorov-Smirnov (KS) distance between two 
consecutive fits of the prior is smaller than \code{tol} AND 
the KS distance between two consecutive fits of the random 
effects variance is smaller than \code{tolrand} OR 2c. The relative difference (Total log marg lik - Previous total log marg lik)/
abs(Total log marg lik) < mliktol. 

Computing time 
is largely determined by the largest element in \code{ntag} and 
\code{maxiter}. We advise to use at least \code{max(ntag)=1000} 
and \code{max(ntag)=2500} when \code{mixtrand=TRUE}. About 
\code{shrinkfixed} and \code{shrinkaddfixed}: if one of the 
fixed parameters is the primary parameter of interest use 
\code{shrinkfixed} for this one, and, when desired, 
\code{shrinkaddfixed} for an another one. 
}

\value{
A list object containing:
  \item{pmlist }{Final estimates of the parameters of the (fitted) priors}
  \item{ksall}{Kolmogorov-Smirnov distances between consecutive estimated priors}
  \item{paraall}{List object containing all consecutive parameter estimates}
  \item{inputpar}{Input parameters used}   
  \item{addfixed}{Parameter values of additional fixed parameters. Passed on for convenience.}
  \item{addrandom}{Parameter values of additional random effects parameters. Passed on for convenience.}
  \item{typelik}{Type of likelihood used}
}

\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.

Van de Wiel MA,  De Menezes RX, Siebring-van Olst E, Van Beusechem VW 
(2013). Analysis of small-sample clinical genomics studies 
using multi-parameter shrinkage: application to high-throughput 
RNA interference screening. BMC Med Genom. 

Van de Wiel MA, Neerincx M, Buffart TE, Sie D, Verheul HMW (2014).
ShrinkBayes: a versatile R-package for analysis of count-based sequencing data in complex study designs. 
BMC Bioinformatics.
} 

\author{
Mark A. van de Wiel
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{ \code{\link[ShrinkBayes]{ShrinkSeq}}

} 

\examples{
\dontrun{
#See Van de Wiel et al. (2013) for more details on the data and the analysis
#To empty current R memory and enlarge memory.limit (Windows) consider, before applying 'ShrinkBayes', using 
#rm(list=ls());gc();
#memory.limit(size = 4000)

#loads HT RNAi screening data; 960 rows (siRNAs), 6 samples
data(HTRNAi)

#Number of cpus to use in (parallel) computations
ncpus2use <- 10

#Covariates treatment and assay, corresponding to the columns of the HTRNAi data. 
treatment <- factor(rep(c("untreated","treated"),3))
assay <- factor(rep(1:3,each=2))

#OFFSET WRT POSITIVE CONTROLS
offsetvalue <- c(0.1703984, -0.6958495,  0.3079694, -0.5582785,  0.2251210, -0.6411269) 

#Formula (as used by INLA)
form = y ~ offset(offsetvalue) + 1 + treatment + assay 

#Simultaneous shrinkage for 'treatment' and 'assay'. In addition, the Gaussian error standard deviation is shrunken by default (see shrinksigma argument). 
shrinksimul <- ShrinkGauss(form=form, dat=HTRNAi,shrinkfixed="treatment",shrinkaddfixed="assay",ncpus=ncpus2use)

#Fit all using the priors resulting from ShrinkGauss
fitg <- FitAllShrink(form,dat=HTRNAi,fams="gaussian",shrinksimul,ncpus=ncpus2use)

#Find nonparametric prior for 'treatment'; note that we explictly allow for a non-symmetric prior (default is symmetric=TRUE)
npprior <- NonParaUpdatePrior(fitall=fitg,modus="fixed",includeP0 = FALSE, shrinkpara="treatment",ncpus=ncpus2use, symmetric=FALSE)

#Update posteriors using the nonparametric prior
nppostshr <- NonParaUpdatePosterior(fitg,npprior,ncpus=ncpus2use)

#Compute local fdrs ( = posterior null (tail-)probabilities under shrinkage prior). 
#Here P(theta <= 0) is computed, assuming {theta<=0} is the null-domain.
lfdr <- SummaryWrap(nppostshr, thr = 0, direction="lesser")

#Compute Bayesian FDRs
BFDRs <- BFDR(lfdr)
}
}



