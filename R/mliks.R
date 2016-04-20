mliks <-function(fital,gauss=TRUE) {
    if(gauss) {unlist(lapply(fital,function(x) {xt <- try(x$mlik[2]); if(mode(xt)=="try-error" | is.null(xt)) {xt <- NA};return(xt)}))
    } else {unlist(lapply(fital,function(x) {xt <- try(x$mlik[1]); if(mode(xt)=="try-error" | is.null(xt)) {xt <- NA};return(xt)}))}
}
