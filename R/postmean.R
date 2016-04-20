postmean <-
function(poster){
if(is.null(poster)) return(NA) else {
return(myinla.expectation(function(x) x,poster))
}
}

