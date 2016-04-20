SpikeSlabPrior <- function(fitall,fitall0){

#fitall=fitzinbt;shrinklc = "MminPcorr"; tol=0.01;maxiter=10;ntotals=c(500,1000);
#modus="fixed"; shrinkpara="group"; ncpus=1;
#shrinklc=FALSE;lcabscoef=c(1,1); vlc=c("fixed","fixed");ndrawtot = 10000;priornew=NULL;p0init=0.7; p0lower = 0.5; pointmass = 0;
#gaussinit0=NULL; gammainit0 = NULL; maxsupport=5; plotdens=FALSE;unimodal=TRUE;logconcave=TRUE;symmetric=FALSE; allow2modes=TRUE;
#fitall <- fitgvague;fitall0 <- fitg0;shrinkpara="group";shrinklc=NULL; ncpus=2;modus="fixed"
#fitall <- ShrinkSeq.fitvague;fitall0=ShrinkSeq.fitzinb0;
fitall <- fitall[[1]]
fitall0 <- fitall0[[1]]

mlik <- mliks(fitall)
mlik0 <- mliks(fitall0)
repNA <- function(x) {if(is.na(x)) return(-10^10) else return(x)}
mlik <- sapply(mlik,repNA)
mlik0 <- sapply(mlik0,repNA)
maxlik <- as.numeric(apply(cbind(mlik,mlik0),1,max))
p0start <- length(which((mlik0-mlik)>0))/length(mlik)

liktot <- function(p0=0.5) {-sum(log(p0*exp(mlik0-maxlik) + (1-p0)*exp(mlik-maxlik)))}
#res1 <- mle(liktot,start=list(p0=p0start),lower=c(0.001),upper=c(0.999),hessian=FALSE,method="L-BFGS-B")
res2 <- optimize(liktot,lower=0,upper=1, maximum = FALSE) #seems more stable than mle
p0est <- res2$minimum
return(p0est)
} #END FUNCTION
