
#setwd("C:\VUData\Maarten\Final");load("MminP_res.Rdata");load("fitzinbTissM.Rdata")
#load("C:\\Synchr\\RPackages\\ShrinkBayes\\Examples\\OutputExamples\\outputExampleScriptsCAGE_1K.Rdata")

plotPoster <- function(id,poster,miny = 10^(-4),xlabel="",ylabel="",plotlabels=NULL,cols = NULL,ltys=NULL,lwds =NULL,xran=NULL,xscatter=NULL,legendxy=NULL,legendbox="n",include0 = TRUE){
#id <- 427;poster <- mixtpostshr;include0 = TRUE;xran<-NULL;miny = 10^(-4);legendbox="n";cols = NULL;ltys=NULL;legendxy=NULL
getpost <- poster[[id]]
postnon0 <- getpost$postbetanon0
if(is.null(plotlabels)) nams <- names(postnon0) else nams <- plotlabels
ell <- length(postnon0)
if(is.null(nams)) nams <- sapply(1:ell,function(i) paste("parameter",i))
rangexy <- function(poster) {
y <- poster[,2]
x <- poster[,1]
sely <- which(y>miny)
xysel <- c(min(x[sely],na.rm=T),max(x[sely],na.rm=T),max(y[sely],na.rm=T))
return(xysel)
}
if(include0){
post0 <- getpost$postbeta0
if(is.null(post0)) post0 <- rep(0,ell)
max0 <- max(post0)
el0 <- length(unique(post0))
} else el0 <- 0
xysels <- sapply(postnon0,rangexy)
xrange <- c(min(xysels[1,]),max(xysels[2,]))
if(!is.null(xran)) xrange <- xran
yrange <- c(0,max(xysels[3,])+0.05*max(xysels[3,]))
if(include0){
yrange[2] <- max(yrange[2], min(1,max0+0.05))
}
if(is.null(cols)) cols <- 1:ell
if(is.null(ltys)) ltys <- rep(1,ell)
if(is.null(lwds)) lwds <- rep(1,ell)

plot(postnon0[[1]],xlim=xrange,ylim=yrange,xlab=xlabel,ylab=ylabel,col=cols[1],lty=ltys[1],lwd=lwds[1],type="l")
if(el0>1) points(c(0,0),c(0,post0[1]),type="l",col=cols[1],lty=ltys[1],lwd=lwds[1])
if(el0==1) points(c(0,0),c(0,post0[1]),type="l",col=cols[1],lty=1,lwd=max(lwds)+2)
if(is.null(xscatter)) xscatter <- (xrange[2]-xrange[1])/200
if(ell>1){
for(i in 2:ell){
xi <- (-1)^i*(i/2)*xscatter
points(postnon0[[i]],xlim=xrange,ylim=yrange,col=cols[i],lty=ltys[i],lwd=lwds[i],type="l")
if(include0 & el0>1) points(c(xi,xi),c(0,post0[i]),type="l",col=cols[i],lty=ltys[i],lwd=lwds[i])
}
}
if(is.null(legendxy)) legendxy <- c(xrange[1],yrange[2] + yrange[2]*0.05)
if(el0>1 | el0==0) legend(legendxy[1],legendxy[2],nams,col=cols,lty=ltys,bty=legendbox) else 
legend(legendxy[1],legendxy[2],c(nams,"pi0All"),col=c(cols,1),lty=c(ltys,1),lwd=c(lwds,max(lwds)+2),bty=legendbox)
}
