fdrcomp <-
function(positives,sig){  #True FDR = FP/N_P. Est FDR = \sum (p0s*I) / N_P
#positives <-  1998;sig <- p0s
sortsig <- sort(sig,index.return=T)
sortind <- sortsig$ix
n <- length(sig)
arr <- rep(0,n)
wh <- which(sortind > positives)
arr[wh] <- 1
FP <- cumsum(arr)
TrueFDR <- FP/1:n
EstFDR <- cumsum(sortsig$x)/(1:n)
return(cbind(TrueFDR,EstFDR))
}

