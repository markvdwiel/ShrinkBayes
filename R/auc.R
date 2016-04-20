auc <-
function(positives,sig,thr=0.2){  #FPR = FP/(N_Neg) = (Pos and H0)/N_H0 = 1-specificity, #TPR = (Pos and H1)/N_Pos = sensitivity
#positives <-  2000;sig <- p0snew
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
whthr <- which(FPR<=thr)
FPRthr <- FPR[whthr]
FPRdif <- FPRthr - c(0,FPRthr[-length(FPRthr)])
TPRthr <- TPR[whthr]
AUC <- TPRthr %*% FPRdif
return(c(AUC=AUC,AUCrel=AUC/(thr^2*.5)))
}

