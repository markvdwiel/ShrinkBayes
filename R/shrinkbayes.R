ShrinkBayesWrap <- function(data,form,paramtotest=NULL,allcontrasts=FALSE,multvscontrol=FALSE,fams=NULL,notfit=NULL,ncpus2use=8,
priorsimple = FALSE,approx0=TRUE,diffthr=0,direction="two-sided",saveposteriors = TRUE, fileposteriors="Posteriors.Rdata", sparse=FALSE,...){
#testContrasts = FALSE, implies K-sample test for factor with more thatn 2 levels
#paramtotest: defaults to first parameter in the model
#notfit: indices of those features for which inference is not desired (e.g. because they did not survive a pre-filter)
#other parameters to feed to ShrinkSeq or ShrinkGauss
#Number of cpus to use in (parallel) computations
#form: right side of formula object only
#group <- factor(c(rep("group1",4),c(rep("group2",4))));data(datsim);data <- datsim;form = ~ 1 + group;paramtotest=NULL;multvscontrol<-FALSE;
#allcontrasts=FALSE;notfit=NULL;ncpus2use=4;approx0=TRUE;saveposteriors<-TRUE;fileposteriors<-"posteriors.RData";diffthr=0;direction="two-sided"
#priorsimple <- TRUE
#...:further arguments to pass on to ShinkSeq or ShrinkGauss
if(sparse) fixed <- c(0,1)
if(allcontrasts & multvscontrol){
    print("Cannot perform all-pairwise and multiple groups versus control comparisons simultaneously")
    cat("Set either \'allcontrasts\' or \'multvscontrol\' to FALSE\n")
    return(NULL)
    }
form <- formula(paste("y",paste(as.character(form)[[1]],as.character(form)[[2]])))

if(is.null(paramtotest)){ 
frmchr <- as.character(form)[[3]]
sp <- strsplit(frmchr,"\\+")
sp <- as.vector(sapply(sp,function(tt){## remove whitespace
    gsub(" ","",tt)
    }))
if(sp[1]=="0" || sp[1]=="1") paramtotest<-sp[2] else paramtotest <- sp[1]
}
print(paste("Performing inference (testing) for parameter(s):",paramtotest))

if(multvscontrol){
print("Performing multiple comparisons with a control group")
print(paste("The control group is:",levels(get(paramtotest))[1]))
cat("If you would like to use a different control group, apply the \'BaselineDef\' function first\n")
}

if(allcontrasts) {
print("Performing inference for all contrasts (using symmetric prior and approximation to null-model)")
lincombvec <- AllComp(paramtotest)
shrinklc <- names(lincombvec)
symmetric <- TRUE
approx0 <- TRUE
} else {
lincombvec <- NULL
shrinklc <- NULL
symmetric <- FALSE
}

Ksam <- FALSE

if(!allcontrasts & !multvscontrol){
if(is.factor(get(paramtotest))){
nlev <- length(levels(get(paramtotest)))
if(nlev > 2){
Ksam <- TRUE
print("Performing K-sample testing")
} else {
print("Performing two-sample test")
}
} else{
print("Performing regression-based test")
}
}

form0 = replacerand0(form,paramtotest) 

datasum <- sum(data[1:5,])
if(is.wholenumber(datasum)) {
counts <- TRUE 
if(is.null(fams)){
print("Data are counts. Zero-inflated negative binomial count model is used")
cat("Use \'fams\' argument to specify a different count model\n")
fams <- "zinb"
}} else {
counts <- FALSE
print("Data are on continuous scale. Gaussian model is used")
fams <- "gaussian"
}

if(Ksam){
excludefornull <- paramtotest
direction <- "equal"
diffthr <- 0
cat("Argument \'direction\' is set to \"equal\" to allow K-sample testing\n")
} else excludefornull <- NULL

print("STARTING INITIAL SHRINKAGE")
pmtinit <- proc.time()
if(counts){
    if(!sparse) shrinksimul <- ShrinkSeq(form=form,dat=data,fams=fams,shrinkfixed=paramtotest, ncpus=ncpus2use,excludefornull=excludefornull,...) else 
      shrinksimul <- ShrinkSeq(form=form,dat=data,fams=fams,shrinkfixed=NULL, ncpus=ncpus2use,excludefornull=excludefornull,fixed=fixed,...)  
    } else {
    if(!sparse) shrinksimul <- ShrinkGauss(form=form, dat=data,shrinkfixed=paramtotest, ncpus=ncpus2use,excludefornull=excludefornull,...) else
      shrinksimul <- ShrinkGauss(form=form,dat=data,fams=fams,shrinkfixed=NULL, ncpus=ncpus2use,excludefornull=excludefornull,fixed=fixed,...)  
    }
time1 <- proc.time()-pmtinit
print("Computing time for shrinkage:")
print(time1)

nr <- nrow(data)
nrfit <- nr - length(notfit)


if(!Ksam){
    if(!is.null(notfit) & (nr > 10000)){
        set.seed(34523)
        rows <- sample(1:nr,size=10000)
        skip<-FALSE
    } else {
        rows <- 1:nr
        skip <- TRUE
        }

print("STARTING INITIAL FIT")
pmt <- proc.time()
    fitg <- FitAllShrink(form,dat=data[rows,],fams=fams,shrinksimul,ncpus=ncpus2use,lincomb=lincombvec)
    
    if(!approx0){
        fitg0 <- FitAllShrink(form0,dat=data[rows,],fams=fams,shrinksimul,ncpus=ncpus2use)
    } else {
        print("Approximating marginal likelihood for null model by Savage-Dickey")
        cat("Set \'approx0 <- FALSE\' when this is not desired.\n")
        fitg0 <- NULL
    }
time1 <- proc.time()-pmt
print("Computing time for initial fit:")
print(time1)  

print("START FITTING MIXTURE PRIOR")
pmt <- proc.time()
if(sparse) prior <- MixtureUpdatePrior(fitall=fitg,fitall0=fitg0, modus="laplace", shrinkpara=paramtotest,shrinklc=shrinklc,ncpus=ncpus2use) else
    if(!priorsimple) prior <- MixtureUpdatePrior(fitall=fitg,fitall0=fitg0, modus="mixt", shrinkpara=paramtotest,shrinklc=shrinklc,ncpus=ncpus2use,symmetric=symmetric) else prior <- MixtureUpdatePrior(fitall=fitg,fitall0=fitg0, modus="gauss", shrinkpara=paramtotest,shrinklc=shrinklc,ncpus=ncpus2use)
time1 <- proc.time()-pmt
print("Computing time for fitting mixture prior:")
print(time1) 

print("START FITTING FOR ALL FEATURES")
pmt <- proc.time()  
    if(!is.null(notfit) & !skip){
        fitg <- FitAllShrink(form,dat=data[-notfit,],fams=fams,shrinksimul,ncpus=ncpus2use,lincomb=lincombvec)
        if(!approx0 | allcontrasts){
            fitg0 <- FitAllShrink(form0,dat=data[-notfit,],fams=fams,shrinksimul,ncpus=ncpus2use)
        } else {
            print("Approximating marginal likelihood for null model by Savage-Dickey")
            cat("Set \'approx0 <- FALSE\' when this is not desired.\n")
            fitg0 <- NULL
        }
        }
    
#save(fitg,prior,fitg0,file="fits.Rdata")

posteriors <- MixtureUpdatePosterior(fitg,prior,fitg0,ncpus=1)
time1 <- proc.time()-pmt
print("Computing time for fitting all features:")
print(time1) 
} else { #K-sample inference
    if(is.null(notfit)) rows <- 1:nr else rows <- (1:nr)[-notfit]
    print("START FITTING FOR ALL FEATURES")
    pmt <- proc.time()  
    fitg <- FitAllShrink(form,dat=data[rows,],fams=fams,shrinksimul,ncpus=ncpus2use,finalprior=TRUE)
    fitg0 <- FitAllShrink(form0,dat=data[rows,],fams=fams,shrinksimul,ncpus=ncpus2use,finalprior=TRUE)
    posteriors <- BFUpdatePosterior(fitg,shrinksimul,fitg0)
    prior <- NULL
    time1 <- proc.time()-pmt
    print("Computing time for fitting all features:")
    print(time1) 
}
print("START COMPUTING SUMMARY STATISTICS, INCL FDR")
pmt <- proc.time() 
if(saveposteriors) save(posteriors,file=fileposteriors)
#check names
res <- SummaryTable(posteriors,BFDRthr=1,diffthr = diffthr,direction=direction,pointmass=0,ndigit=3,ncpus=1)
if(!is.null(notfit)) res$index <- (1:nr)[-notfit][res$index]
cn <- colnames(res)
whBFDR <- grep("BFDR",cn)
nsigsFDR01 <- sapply(whBFDR,function(wh) length(which(res[,wh]<=0.1)))
names(nsigsFDR01) <- cn[whBFDR]
print("Number of significant results at (B)FDR <= 0.1:")
print(nsigsFDR01)
time1 <- proc.time()-pmt
print("Computing time for summary statistics:")
print(time1) 
time2 <- proc.time()-pmtinit
print("Total computing time:")
print(time2) 
return(list(FDRs = res, nsigsFDR01=nsigsFDR01, shrink = shrinksimul, prior = prior))
}
