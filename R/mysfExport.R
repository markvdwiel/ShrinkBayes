mysfExport <-
function(low=10,forceexport=NULL,forceexclude=NULL){
#low: upper threshold on size of objects (in Mb) to be exported by default
#forceexport: character vector with names of (large) objects that should always be exported no matter how large they are
objsize <- sapply(objects(envir = globalenv()),function(obj) {gobj <- get(obj);object.size(gobj)})
whlarge <- which(objsize > low*10^6)
if(length(whlarge) > 0){
whexp <- names(objsize)[whlarge]
notexp <- union(setdiff(whexp,forceexport),forceexclude)
sfExportAll(except=notexp)
} else {
sfExportAll()
}
}

