samplepost <-
function(ex,reps){
sampost <- try(myinla.rmarginal(reps,ex),silent=TRUE)
if(class(sampost)=="try-error") {
    yprobs <- ex[,2]/sum(ex[,2])
    sampost <- try(sample(ex[,1],size=reps,prob=yprobs,replace=TRUE), silent=TRUE)
    if(class(sampost)=="try-error") sampost <- NA   
}
return(sampost)
}

