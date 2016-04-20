mixture2 <-
function(x,thrfit = 8){
    p0 <- sum(abs(x)==Inf)/length(x)
    x <- x[abs(x) <= thrfit]
    out <- list(mixp = c(p0,1-p0), meanprec = c(mean(x),1/(sd(x))^2))
    return(out)
}

