roc <-
function(positives,sig){  #FPR = FP/(N_Neg) = (Pos and H0)/N_H0 = 1-specificity, #TPR = (Pos and H1)/N_Pos = sensitivity
#positives <-  1998;sig <- p0s
sortsig <- sort(sig,index.return=T)
sortind <- sortsig$ix
n <- length(sig)
arr <- rep(0,n)
wh <- which(sortind > positives)
arr[wh] <- 1
FP <- cumsum(arr)
FPR <- FP/(n-positives)
TP <- (1:n)-FP
TPR <- TP/positives
return(cbind(FPR,TPR))
}

