lfdr <-
function(poster, thr,direction){
if(is.null(poster)) return(NA) else {
marg <- try(myinla.pmarginal(thr,poster))
if(class(marg)=="try-error") return(NA) else
    {
    if(direction=="greater") toret <- 1-marg else toret <- marg
    return(toret)
    }
}
}
