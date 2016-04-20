FastScreenP <-
function(pvalues,method="padjust",threshold=ifelse(method=="pvalue",0.1,0.5),adjmethod="BH"){
#method=padjust or method==pvalue
if(method=="padjust") pvalues <- p.adjust(pvalues,method=adjmethod)
wh <- which(pvalues <= threshold)
return(wh)
}

