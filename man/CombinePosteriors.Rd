\name{CombinePosteriors}
\alias{CombinePosteriors}
\title{
Combines posteriors from two models
}
\description{
This function combines posteriors from two models for a given parameter and, possibly, for related linear combinations. The mixture prior probabilities 
should be present in \code{shrinksimul}.
}
\usage{
CombinePosteriors(fullpost1, fullpost2, shrinksimul, para, includelincomb = TRUE, margcomb = c("marginals.fixed"), ncpus = 2, ndigits = 3, updateby = 1000)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fullpost1}{A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}}
}
  \item{fullpost2}{
A 2-component list object resulting from \code{\link[ShrinkBayes]{FitAllShrink}}
}
  \item{shrinksimul}{
A list object resulting from \code{\link[ShrinkBayes]{ShrinkSeq}} or \code{\link[ShrinkBayes]{ShrinkGauss}}
}
  \item{para}{Character string containing the name of the parameter for which combined posteriors are desired. 
}
  \item{includelincomb}{
Boolean indicating whether combined posteriors for linear combinations are desired.
}
  \item{margcomb}{Character string. Either "marginals.fixed", "marginals.random", "marginals.hyper" or "marginals.internal.hyper".
}
  \item{ncpus}{Integer. The number of cpus to use for parallel computations.
}
  \item{ndigits}{Integer. Numerical precision in digits for the output.
}
  \item{updateby}{Integer. Perform the computation in blocks of size \code{updateby}. In particular useful when fullpost1 and fullpost2 are large objects.
}
}
\details{
The parameter names of hyper-parameters like zero-inflation, random effects precision and size of the negative binomial are automatically generated 
by \code{INLA}. The names may be somewhat incovenient. Please check these, e.g. by using \code{names(fullpost1[[1]][[1]]$marginals.hyper)}.
}
\value{
Two-component list object.
%%  If it is a LIST, use
  \item{res}{A list that contains the combined posteriors in numerical format  }
  \item{priors}{A list passing through all information on priors present in \code{fullpost1}}
%% ...
}

\references{
Van de Wiel MA, Leday GGR, Pardo L, Rue H, Van der Vaart AW, Van Wieringen WN (2012). 
Bayesian analysis of RNA sequencing data by estimating multiple shrinkage priors. Biostatistics.
}

\author{
Mark A. van de Wiel
}

\seealso{
\code{\link[ShrinkBayes]{FitAllShrink}}, \code{\link[ShrinkBayes]{ShrinkSeq}}, \code{\link[ShrinkBayes]{ShrinkGauss}}, 
}

\examples{
\dontrun{
#See Van de Wiel et al. (2012) for more details on the data and the analysis
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
shrinksimul <- ShrinkSeq(form=form, dat=CAGEdata,shrinkfixed="groupfac",shrinkrandom="pers",mixtdisp=TRUE,ncpus=ncpus2use, maxiter=3)

#The mass on the point mass of the mixture prior for overdispersion; 
#if close to zero, then only fits on (zi-)-negative binomial are needed; otherwise (zi-)Poisson is also needed
shrinksimul$pmlist$mixp

#Fit all using the priors resulting from ShrinkSeq and zi-poisson likelihood; also find posteriors for the linear combinations defined above
fitzip <- FitAllShrink(form,dat=CAGEdata,fams="zip",shrinksimul,ncpus=ncpus2use,lincomb=lincombvec)

#Fit all using the priors resulting from ShrinkSeq and zi-nb likelihood; also find posteriors for the linear combinations defined above
fitzinb <- FitAllShrink(form,dat=CAGEdata,fams="zinb",shrinksimul,ncpus=ncpus2use,lincomb=lincombvec)

#Combine posteriors from fitzip and fitzinb; linear combinations are included by default
cp <- CombinePosteriors(fitzip,fitzinb,shrinksimul,para="groupfac",ncpus=ncpus2use) 

#Find nonparametric prior for all group-wise differences. These include the linear combinations. 
#Please check lcfac argument (defaults to 2): should be sum of squares of weights in the lincombvec (= 2 for pairwise contrasts)
#INCREASE maxiter FOR OBTAINING FINAL RESULTS (DEFAULT = 15)
npprior <- NonParaUpdatePrior(fitall=cp,modus="fixed", shrinkpara="groupfac", shrinklc=TRUE,ncpus=ncpus2use, maxiter = 3)

#Update posteriors using the nonparametric prior
nppostshr <- NonParaUpdatePosterior(cp,npprior,ncpus=ncpus2use)

#Compute local fdrs ( = posterior null (tail-)probabilities under shrinkage prior)
#Here P(theta <= 0.1) is computed, assuming {theta<=0.1} is the null-domain.
lfdr <- SummaryWrap(nppostshr, thr = 0.1)

#Compute Bayesian FDRs; all BFDRs for each comparison [matrix format]
BFDRs <- BFDR(lfdr)

#Compute Bayesian FDRs; BFDR for multiple comparisons
BFDRmult <- BFDR(lfdr,multcomp=TRUE)
}
}




