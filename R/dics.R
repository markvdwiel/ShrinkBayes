dics <-
function(fital) {unlist(lapply(fital,function(x) if(!is.null(x$dic[1])) x$dic[1] else NA))}

