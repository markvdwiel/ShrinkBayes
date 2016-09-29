### R code from vignette source 'ShrinkBayes.Rnw'

###################################################
### code chunk number 1: ShrinkBayes.Rnw:129-132
###################################################
library(ShrinkBayes)
data(mirseqnorm)
head(mirseqnorm[,1:10])


###################################################
### code chunk number 2: ShrinkBayes.Rnw:137-139
###################################################
data(designmirseq)
head(designmirseq)


###################################################
### code chunk number 3: ShrinkBayes.Rnw:144-149
###################################################
PM <- designmirseq$PM
indiv <- designmirseq$indiv
timepos <- designmirseq$timepos
chemo <- designmirseq$chemo
organ <- designmirseq$organ


###################################################
### code chunk number 4: ShrinkBayes.Rnw:154-155
###################################################
form = ~ 1 + PM + timepos + chemo + organ + f(indiv)


###################################################
### code chunk number 5: ShrinkBayes.Rnw:160-161 (eval = FALSE)
###################################################
## SBmir <- ShrinkBayesWrap(mirseqnorm,form)  


###################################################
### code chunk number 6: ShrinkBayes.Rnw:171-173 (eval = FALSE)
###################################################
## SBmirsmall <- ShrinkBayesWrap(mirseqnorm[1:100,],form,
## ntag=c(25,50),maxiter=2,priorsimple=TRUE, approx0=TRUE)  


###################################################
### code chunk number 7: ShrinkBayes.Rnw:181-183 (eval = FALSE)
###################################################
## SBmirshrink <- ShrinkBayesWrap(mirseqnorm,form,
## shrinkaddfixed=c("organ","chemo","timepos"))  


###################################################
### code chunk number 8: ShrinkBayes.Rnw:189-190 (eval = FALSE)
###################################################
## SBmirorgan <- ShrinkBayesWrap(mirseqnorm,form,paramtotest="organ")  


###################################################
### code chunk number 9: ShrinkBayes.Rnw:214-217
###################################################
library(ShrinkBayes)
data(datsim)
head(datsim[,1:5])


###################################################
### code chunk number 10: ShrinkBayes.Rnw:222-223
###################################################
ncpus2use <- 10


###################################################
### code chunk number 11: ShrinkBayes.Rnw:231-232
###################################################
group <- factor(c(rep("group1",4),c(rep("group2",4))))


###################################################
### code chunk number 12: ShrinkBayes.Rnw:237-238
###################################################
form = y ~  1 + group 


###################################################
### code chunk number 13: ShrinkBayes.Rnw:244-246 (eval = FALSE)
###################################################
## shrinksimul <- ShrinkGauss(form=form, dat=datsim,shrinkfixed="group", 
## ncpus=ncpus2use) 


###################################################
### code chunk number 14: ShrinkBayes.Rnw:257-258
###################################################
load("C:\\Synchr\\RPackages\\ShrinkBayes\\Examples\\OutputExamples\\outputToyExample.Rdata")


###################################################
### code chunk number 15: ShrinkBayes.Rnw:262-264 (eval = FALSE)
###################################################
## fitg <- FitAllShrink(form,dat=datsim,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)


###################################################
### code chunk number 16: ShrinkBayes.Rnw:271-272
###################################################
form0 = y ~  1 


###################################################
### code chunk number 17: ShrinkBayes.Rnw:277-279 (eval = FALSE)
###################################################
## fitg0 <- FitAllShrink(form0,dat=datsim,fams="gaussian",shrinksimul,
## ncpus=ncpus2use) 


###################################################
### code chunk number 18: ShrinkBayes.Rnw:294-297 (eval = FALSE)
###################################################
## mixtprior2gauss <- 
## MixtureUpdatePrior(fitall=fitg,fitall0=fitg0, modus="mixt", 
## shrinkpara="group",ncpus=ncpus2use) 


###################################################
### code chunk number 19: ShrinkBayes.Rnw:300-302
###################################################
bestfinal <- mixtprior2gauss$allpara[1,]
bestfinal


###################################################
### code chunk number 20: ShrinkBayes.Rnw:305-307 (eval = FALSE)
###################################################
## mixtpostshr <- MixtureUpdatePosterior(fitg,mixtprior2gauss,fitg0,
## ncpus=ncpus2use)


###################################################
### code chunk number 21: ShrinkBayes.Rnw:311-313 (eval = FALSE)
###################################################
## lfdrless <- SummaryWrap(mixtpostshr, thr = 0, direction="lesser")
## lfdrgreat <- SummaryWrap(mixtpostshr, thr = 0, direction="greater")


###################################################
### code chunk number 22: ShrinkBayes.Rnw:324-325 (eval = FALSE)
###################################################
## BFDRs <- BFDR(lfdrless,lfdrgreat)


###################################################
### code chunk number 23: ShrinkBayes.Rnw:330-331
###################################################
plot(BFDRs)


###################################################
### code chunk number 24: ShrinkBayes.Rnw:347-359
###################################################
fdrcomp <- function(negatives,sig){  #True FDR = FP/N_P. Est FDR = \sum (p0s*I) / N_P
sortsig <- sort(sig,index.return=T)
sortind <- sortsig$ix
n <- length(sig)
arr <- rep(0,n)
wh <- which(sapply(sortind,is.element,set=negatives))
arr[wh] <- 1
FP <- cumsum(arr)
TrueFDR <- FP/1:n
EstFDR <- sortsig$x
return(cbind(TrueFDR,EstFDR))
}


###################################################
### code chunk number 25: ShrinkBayes.Rnw:365-368
###################################################
res <- fdrcomp(1:1000,BFDRs)
plot(res,type="l")
abline(a=0,b=1,col="red")


###################################################
### code chunk number 26: ShrinkBayes.Rnw:374-376 (eval = FALSE)
###################################################
## prior1Normalp0 <- MixtureUpdatePrior(fitall=fitg,fitall=fitg0, 
## modus="gauss", shrinkpara="group",ncpus=ncpus2use)


###################################################
### code chunk number 27: ShrinkBayes.Rnw:386-388
###################################################
bestfinal <- prior1Normalp0$allpara[1,]
bestfinal


###################################################
### code chunk number 28: ShrinkBayes.Rnw:395-400 (eval = FALSE)
###################################################
## post1Normalp0shr <- MixtureUpdatePosterior(fitg,prior1Normalp0,fitg0,
## ncpus=ncpus2use)
## lfdrless2 <- SummaryWrap(post1Normalp0shr, thr = 0, direction="lesser")
## lfdrgreat2 <- SummaryWrap(post1Normalp0shr, thr = 0, direction="greater")
## BFDRs2 <- BFDR(lfdrless2,lfdrgreat2)


###################################################
### code chunk number 29: ShrinkBayes.Rnw:406-408
###################################################
findef <- par("fin")
par(fin=c(4,4))


###################################################
### code chunk number 30: ShrinkBayes.Rnw:411-412
###################################################
plot(BFDRs,BFDRs2,type="l")


###################################################
### code chunk number 31: ShrinkBayes.Rnw:420-421
###################################################
par(fin=findef)


###################################################
### code chunk number 32: ShrinkBayes.Rnw:425-428
###################################################
library(ShrinkBayes)
data(HTRNAi)
head(HTRNAi)


###################################################
### code chunk number 33: ShrinkBayes.Rnw:433-434
###################################################
ncpus2use <- 10


###################################################
### code chunk number 34: ShrinkBayes.Rnw:439-441
###################################################
treatment <- factor(rep(c("untreated","treated"),3))
assay <- factor(rep(1:3,each=2))


###################################################
### code chunk number 35: ShrinkBayes.Rnw:446-448
###################################################
offsetvalue <- c(0.1703984, -0.6958495,  0.3079694, -0.5582785,  
0.2251210, -0.6411269)


###################################################
### code chunk number 36: ShrinkBayes.Rnw:457-458
###################################################
form = y ~ offset(offsetvalue) + 1 + treatment + assay


###################################################
### code chunk number 37: ShrinkBayes.Rnw:463-465 (eval = FALSE)
###################################################
## shrinksimul <- ShrinkGauss(form=form, dat=HTRNAi,shrinkfixed="treatment",
## shrinkaddfixed="assay", fixedmeanzero = FALSE, ncpus=ncpus2use)


###################################################
### code chunk number 38: ShrinkBayes.Rnw:494-495
###################################################
load("C:\\Synchr\\RPackages\\ShrinkBayes\\Examples\\OutputExamples\\outputExampleScriptsgauss.Rdata")


###################################################
### code chunk number 39: ShrinkBayes.Rnw:499-500
###################################################
shrinksimul$pmlist


###################################################
### code chunk number 40: ShrinkBayes.Rnw:506-507
###################################################
round(shrinksimul$paraall[,1:6],3)


###################################################
### code chunk number 41: ShrinkBayes.Rnw:510-511
###################################################
round(shrinksimul$paraall[,-(1:6)],3)


###################################################
### code chunk number 42: ShrinkBayes.Rnw:517-519 (eval = FALSE)
###################################################
## fitg <- FitAllShrink(form,dat=HTRNAi,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)


###################################################
### code chunk number 43: ShrinkBayes.Rnw:542-543
###################################################
fit1 <- fitg$res[[1]]


###################################################
### code chunk number 44: ShrinkBayes.Rnw:550-551
###################################################
fit1$summary.fixed


###################################################
### code chunk number 45: ShrinkBayes.Rnw:555-557
###################################################
marginal <- fit1$marginals.fixed$treatmentuntreated
plot(marginal)


###################################################
### code chunk number 46: ShrinkBayes.Rnw:564-565
###################################################
fit1$summary.hyper


###################################################
### code chunk number 47: ShrinkBayes.Rnw:570-571
###################################################
fit1$mlik


###################################################
### code chunk number 48: ShrinkBayes.Rnw:575-578 (eval = FALSE)
###################################################
## npprior <- NonParaUpdatePrior(fitall=fitg,modus="fixed", 
## shrinkpara="treatment",ncpus=ncpus2use, includeP0=FALSE, 
## allow2modes=FALSE, symmetric=FALSE, logconcave=TRUE)


###################################################
### code chunk number 49: ShrinkBayes.Rnw:591-595
###################################################
plot(npprior$priornew,type="l")
supp <- npprior$priornew[,1]
points(supp,dnorm(supp,mean=shrinksimul$pmlist$mufixed, 
sd  = 1/sqrt(shrinksimul$pmlist$precfixed)),type="l",col="red",lty=2)


###################################################
### code chunk number 50: ShrinkBayes.Rnw:601-602 (eval = FALSE)
###################################################
## nppostshr <- NonParaUpdatePosterior(fitg,npprior,ncpus=ncpus2use)


###################################################
### code chunk number 51: ShrinkBayes.Rnw:606-607 (eval = FALSE)
###################################################
## lfdr <- SummaryWrap(nppostshr, thr = 0, direction="lesser")


###################################################
### code chunk number 52: ShrinkBayes.Rnw:613-614 (eval = FALSE)
###################################################
## BFDRs <- BFDR(lfdr) 


###################################################
### code chunk number 53: ShrinkBayes.Rnw:619-629
###################################################
whsig <- which(BFDRs <= 0.1)
whsig
BFDRs[whsig]
layout(matrix(1:3,nrow=1))
plot(nppostshr[[whsig[1]]][[1]][[1]],xlim=c(-0.5,0.7),type="l")
abline(v=0,lty=2)
plot(nppostshr[[whsig[2]]][[1]][[1]],xlim=c(-0.5,0.7),type="l")
abline(v=0,lty=2)
plot(nppostshr[[whsig[3]]][[1]][[1]],xlim=c(-0.5,0.7),type="l")
abline(v=0,lty=2)


###################################################
### code chunk number 54: ShrinkBayes.Rnw:637-638
###################################################
layout(matrix(c(1),nrow=1))


###################################################
### code chunk number 55: ShrinkBayes.Rnw:646-651
###################################################
library(ShrinkBayes)
data(CAGEdata10000)
CAGEdata <- CAGEdata10000
CAGEdata <- CAGEdata[1:1000,]
CAGEdata[1:2,]


###################################################
### code chunk number 56: ShrinkBayes.Rnw:655-657
###################################################
data(design_brain)
design_brain


###################################################
### code chunk number 57: ShrinkBayes.Rnw:661-664
###################################################
pers <- design_brain$pers  #persons
batch <-design_brain$batch   #batch
groupfac <- design_brain$groupfac #group (= brain region)


###################################################
### code chunk number 58: ShrinkBayes.Rnw:668-669
###################################################
ncpus2use <- 10


###################################################
### code chunk number 59: ShrinkBayes.Rnw:674-675
###################################################
groupfac <- BaselineDef("groupfac",baselinegroup="1")


###################################################
### code chunk number 60: ShrinkBayes.Rnw:685-686
###################################################
lincombvec <- AllComp("groupfac")


###################################################
### code chunk number 61: ShrinkBayes.Rnw:689-690
###################################################
form = y ~ 1 + groupfac + batch + f(pers,model="iid") 


###################################################
### code chunk number 62: ShrinkBayes.Rnw:696-698 (eval = FALSE)
###################################################
## shrinksimul <- ShrinkSeq(form=form, dat=CAGEdata,shrinkfixed="groupfac",
## shrinkrandom="pers",mixtdisp=TRUE,ncpus=ncpus2use)


###################################################
### code chunk number 63: ShrinkBayes.Rnw:712-713
###################################################
load("C:\\Synchr\\RPackages\\ShrinkBayes\\Examples\\OutputExamples\\outputExampleScriptsCAGE_1K.Rdata")


###################################################
### code chunk number 64: ShrinkBayes.Rnw:717-718
###################################################
shrinksimul$pmlist$mixp


###################################################
### code chunk number 65: ShrinkBayes.Rnw:725-727 (eval = FALSE)
###################################################
## fitzip <- FitAllShrink(form,dat=CAGEdata,fams="zip",shrinksimul,
## ncpus=ncpus2use,lincomb=lincombvec)


###################################################
### code chunk number 66: ShrinkBayes.Rnw:733-735 (eval = FALSE)
###################################################
## fitzinb <- FitAllShrink(form,dat=CAGEdata,fams="zinb",shrinksimul,
## ncpus=ncpus2use,lincomb=lincombvec)


###################################################
### code chunk number 67: ShrinkBayes.Rnw:740-742 (eval = FALSE)
###################################################
## cp <- CombinePosteriors(fitzip,fitzinb,shrinksimul,para="groupfac",
## ncpus=ncpus2use) 


###################################################
### code chunk number 68: ShrinkBayes.Rnw:749-752 (eval = FALSE)
###################################################
## npprior <- NonParaUpdatePrior(fitall=cp,modus="fixed", shrinkpara="groupfac", 
## shrinklc=names(lincombvec),lincombs=lincombvec,ncpus=ncpus2use, maxiter=3, 
## includeP0 = FALSE, symmetric = TRUE, allow2modes=FALSE)


###################################################
### code chunk number 69: ShrinkBayes.Rnw:772-774
###################################################
theprior <- npprior$priornew
plot(theprior,type="l",xlim=c(-2.5,2.5))


###################################################
### code chunk number 70: ShrinkBayes.Rnw:777-783
###################################################
quantiles <- inla.qmarginal(p=c(0.05,0.25,0.5,0.75,0.95), theprior)
quantiles
expect <- inla.emarginal(function(x) x, theprior)
expect
sd <- sqrt(inla.emarginal(function(x) x^2, theprior) - expect^2)
sd


###################################################
### code chunk number 71: ShrinkBayes.Rnw:788-789 (eval = FALSE)
###################################################
## nppostshr <- NonParaUpdatePosterior(cp,npprior,ncpus=ncpus2use)


###################################################
### code chunk number 72: ShrinkBayes.Rnw:795-796 (eval = FALSE)
###################################################
## lfdr <- SummaryWrap(nppostshr, thr = log(1.5))


###################################################
### code chunk number 73: ShrinkBayes.Rnw:807-808 (eval = FALSE)
###################################################
## BFDRs <- BFDR(lfdr)


###################################################
### code chunk number 74: ShrinkBayes.Rnw:810-811
###################################################
head(BFDRs)


###################################################
### code chunk number 75: ShrinkBayes.Rnw:816-817 (eval = FALSE)
###################################################
## BFDRmult <- BFDR(lfdr,multcomp=TRUE)


###################################################
### code chunk number 76: ShrinkBayes.Rnw:825-827
###################################################
wh <- which(BFDRmult <= 0.1)
length(wh)


###################################################
### code chunk number 77: ShrinkBayes.Rnw:830-835
###################################################
whcomp <- which(BFDRs[wh,]<= 0.1,arr.ind=TRUE)
plot(cbind(whcomp[,2],wh[whcomp[,1]]),type="p",xlab="Comparison",
ylab="Feature index",xaxt="n")
axis(1,at=1:10,labels=c("1-2","1-3","1-4","1-5","2-3","2-4","3-4",
"2-5","3-5","4-5"))


###################################################
### code chunk number 78: ShrinkBayes.Rnw:1024-1030 (eval = FALSE)
###################################################
## #1500 rows (siRNAs), 8 samples, pi0 = 2/3
## datsim1 <- matrix(rnorm(8000,mean=0,sd=0.5),nrow=1000)
## meanvec <- matrix(rep(rnorm(500,0,1),4),nrow=500)
## datsim2 <- cbind(matrix(rnorm(2000,mean=0,sd=0.5),nrow=500),
## matrix(rnorm(2000,mean=0,sd=0.5),nrow=500) + meanvec)
## datsim <- rbind(datsim1,datsim2)


###################################################
### code chunk number 79: ShrinkBayes.Rnw:1034-1053 (eval = FALSE)
###################################################
## library(ShrinkBayes)
## data(datsim)
## ncpus2use <- 10
## group <- factor(c(rep("group1",4),c(rep("group2",4))))
## form = y ~  1 + group
## form0 = y ~ 1 
## shrinksimul <- ShrinkGauss(form=form, dat=datsim,shrinkfixed="group", 
## ncpus=ncpus2use) 
## fitg <- FitAllShrink(form,dat=datsim,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## fitg0 <- FitAllShrink(form0,dat=datsim,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## mixtprior2gauss <- MixtureUpdatePrior(fitall=fitg,fitall0=fitg0, 
## modus="mixt", shrinkpara="group",ncpus=ncpus2use)
## mixtpostshr <- MixtureUpdatePosterior(fitg,mixtprior2gauss,fitg0,
## ncpus=ncpus2use)
## lfdrless <- SummaryWrap(mixtpostshr, thr = 0, direction="lesser")
## lfdrgreat <- SummaryWrap(mixtpostshr, thr = 0, direction="greater")
## BFDRs <- BFDR(lfdrless,lfdrgreat)


###################################################
### code chunk number 80: ShrinkBayes.Rnw:1057-1075 (eval = FALSE)
###################################################
## library(ShrinkBayes)
## data(HTRNAi)
## ncpus2use <- 10
## treatment <- factor(rep(c("untreated","treated"),3))
## assay <- factor(rep(1:3,each=2))
## offsetvalue <- c(0.1703984, -0.6958495,  0.3079694, -0.5582785,  
## 0.2251210, -0.6411269)
## form = y ~ offset(offsetvalue) + 1 + treatment + assay
## shrinksimul <- ShrinkGauss(form=form, dat=HTRNAi,shrinkfixed="treatment",
## shrinkaddfixed="assay", fixedmeanzero = FALSE, ncpus=ncpus2use)
## fitg <- FitAllShrink(form,dat=HTRNAi,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## npprior <- NonParaUpdatePrior(fitall=fitg,modus="fixed", 
## shrinkpara="treatment", ncpus=ncpus2use, includeP0 = FALSE, 
## logconcave=TRUE, allow2modes=FALSE)
## nppostshr <- NonParaUpdatePosterior(fitg,npprior,ncpus=ncpus2use)
## lfdr <- SummaryWrap(nppostshr, thr = 0, direction="lesser")
## BFDRs <- BFDR(lfdr) 


###################################################
### code chunk number 81: ShrinkBayes.Rnw:1079-1103 (eval = FALSE)
###################################################
## data(CAGEdata10000)
## CAGEdata <- CAGEdata10000[1:1000,]
## data(design_brain)
## pers <- design_brain$pers ; batch <-design_brain$batch; 
## groupfac <- design_brain$groupfac 
## ncpus2use <- 10
## groupfac <- BaselineDef("groupfac",baselinegroup="1")
## lincombvec <- AllComp("groupfac")
## form = y ~ 1 + groupfac + batch + f(pers,model="iid") 
## shrinksimul <- ShrinkSeq(form=form, dat=CAGEdata,shrinkfixed="groupfac",
## shrinkrandom="pers",mixtdisp=TRUE,ncpus=ncpus2use)
## fitzip <- FitAllShrink(form,dat=CAGEdata,fams="zip",shrinksimul,
## ncpus=ncpus2use,lincomb=lincombvec)
## fitzinb <- FitAllShrink(form,dat=CAGEdata,fams="zinb",shrinksimul,
## ncpus=ncpus2use,lincomb=lincombvec)
## cp <- CombinePosteriors(fitzip,fitzinb,shrinksimul,para="groupfac",
## ncpus=ncpus2use) 
## npprior <- NonParaUpdatePrior(fitall=cp,modus="fixed", 
## shrinkpara="groupfac", shrinklc=TRUE,ncpus=ncpus2use, maxiter=3, 
## includeP0 = FALSE, symmetric=TRUE, logconcave=FALSE, allow2modes=FALSE)
## nppostshr <- NonParaUpdatePosterior(cp,npprior,ncpus=ncpus2use)
## lfdr <- SummaryWrap(nppostshr, thr = log(1.5))
## BFDRs <- BFDR(lfdr)
## BFDRmult <- BFDR(lfdr,multcomp=TRUE)


###################################################
### code chunk number 82: ShrinkBayes.Rnw:1107-1113 (eval = FALSE)
###################################################
## datsim1 <- matrix(rnorm(400000*8,mean=0,sd=0.5),nrow=400000)
## meanvec <- matrix(rep(rnorm(200000,0,1),4),nrow=200000)
## datsim2 <- cbind(matrix(rnorm(200000*4,mean=0,sd=0.5),nrow=200000),
## matrix(rnorm(200000*4,mean=0,sd=0.5),nrow=200000) + meanvec)
## datsimlarge <- rbind(datsim1,datsim2)
## save(datsimlarge,file="datsimlarge.Rdata")


###################################################
### code chunk number 83: ShrinkBayes.Rnw:1118-1120 (eval = FALSE)
###################################################
## tuningsize <- 10000 
## batchsize <- 50000 


###################################################
### code chunk number 84: ShrinkBayes.Rnw:1133-1137 (eval = FALSE)
###################################################
## whrows <- sample(1:nrow(datsimlarge),tuningsize) 
## datsimtune <- datsimlarge[whrows,] 
## save(datsimtune,file="datsimtune.Rdata") 
## rm(datsimlarge);gc() 


###################################################
### code chunk number 85: ShrinkBayes.Rnw:1142-1147 (eval = FALSE)
###################################################
## ncpus2use <- 6
## group <- factor(c(rep("group1",4),c(rep("group2",4))))
## form = y ~  1 + group 
## form0 = y ~  1
## library(ShrinkBayes)


###################################################
### code chunk number 86: ShrinkBayes.Rnw:1151-1160 (eval = FALSE)
###################################################
## shrinksimul <- ShrinkGauss(form=form, dat=datsimtune,shrinkfixed="group", 
## ncpus=ncpus2use) 
## fitg <- FitAllShrink(form,dat=datsimtune,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## fitg0 <- FitAllShrink(form0,dat=datsimtune,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## mixtprior<- MixtureUpdatePrior(fitall=fitg, fitall0=fitg0,modus="gauss", 
## shrinkpara="group",ncpus=ncpus2use)
## save(shrinksimul,mixtprior,file="priorstuned.Rdata")


###################################################
### code chunk number 87: ShrinkBayes.Rnw:1166-1167 (eval = FALSE)
###################################################
## saveFits <- FALSE


###################################################
### code chunk number 88: ShrinkBayes.Rnw:1172-1173 (eval = FALSE)
###################################################
## savePosteriors <- TRUE


###################################################
### code chunk number 89: ShrinkBayes.Rnw:1180-1182 (eval = FALSE)
###################################################
## load("priorstuned.Rdata")
## load("datsimlarge.Rdata")


###################################################
### code chunk number 90: ShrinkBayes.Rnw:1186-1190 (eval = FALSE)
###################################################
## nr <- nrow(datsimlarge)
## nloop <- ceiling(nr/batchsize)
## lfdrless <- c()
## lfdrgreat <- c()


###################################################
### code chunk number 91: ShrinkBayes.Rnw:1194-1220 (eval = FALSE)
###################################################
## for(k in 1:nloop){
## if(k>1) load("datsimlarge.Rdata")
## rangek <- ((k-1)*batchsize+1):min(nr,k*batchsize)
## datk <- datsimlarge[rangek,]
## print(paste("Computing posteriors for features",rangek[1],"to",
## rangek[length(rangek)])) 
## rm(datsimlarge);gc()
## fitgk <- FitAllShrink(form,dat=datk,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## fitg0k <- FitAllShrink(form0,dat=datk,fams="gaussian",shrinksimul,
## ncpus=ncpus2use)
## if(saveFits) {
##     save(fitgk,file=paste("fitg_batch_",k,".Rdata",sep=""))
##     save(fitg0k,file=paste("fitg0_batch_",k,".Rdata",sep=""))
##     }
## mixtpostshrk <- MixtureUpdatePosterior(fitgk,mixtprior,fitg0k,
## ncpus=ncpus2use)
## if(savePosteriors) save(mixtpostshrk,
## file=paste("mixtpostshr_batch_",k,".Rdata",sep=""))
## lfdrlessk <- SummaryWrap(mixtpostshrk, thr = 0, direction="lesser")
## lfdrgreatk <- SummaryWrap(mixtpostshrk, thr = 0, direction="greater")
## lfdrless <- rbind(lfdrless,lfdrlessk)
## lfdrgreat <- rbind(lfdrgreat,lfdrgreatk)
## save(lfdrless,lfdrgreat,file="lfdrs.Rdata")
## rm(fitgk,mixtpostshr);gc()
## }


###################################################
### code chunk number 92: ShrinkBayes.Rnw:1227-1228 (eval = FALSE)
###################################################
## BFDRs <- BFDR(lfdrless,lfdrgreat)


