#randreplace <-
#function(form,shrinkrandom,initrandomprec){
#frmchr <- as.character(form)[[3]]
#sp <- strsplit(frmchr,"\\+")
#sp <- as.vector(sapply(sp,function(tt){## remove whitespace
#    gsub(" ","",tt)
#    }))
#wh <- grep(shrinkrandom,sp)
#if(is.numeric(wh)) {
#sp[wh] <- paste("f(",shrinkrandom,", model =\"iid\", prior=\"loggamma\", param=c(",initrandomprec[1],",",initrandomprec[2],"))",sep="")
#fnew <- sp[1]; if(length(sp)>1) {for(i in 2:length(sp)) fnew <- paste(fnew,sp[i],sep="+")}
#formret <- formula(paste("y ~ ",fnew,sep=""))
#return(formret)
#} else {
#    return(form)
#}
#}

randreplace <-
function(form,shrinkrandom,initrandomprec){
#form <- form2b; shrinkrandom="pers";initrandomprec <- c(10,2)
frmchr <- as.character(form)[[3]]
sp <- strsplit(frmchr,"\\+")
sp <- as.vector(sapply(sp,function(tt){## remove whitespace
    gsub(" ","",tt)
    }))
sr1 <- paste("\\(",shrinkrandom,",",sep="")  
sr2 <- paste("\\(",shrinkrandom,"\\)",sep="")   
wh <- union(grep(sr1,sp),grep(sr2,sp))
if(is.numeric(wh)) {
chr2 <- sp[wh]
sp2 <- strsplit(chr2,",")
sp2 <- as.vector(sapply(sp2,function(tt){## remove whitespace
    gsub(" ","",tt)
    }))
wh2 <- grep("param=",sp2)
if(length(wh2)!=0) {
l2 <- length(grep("\\(",sp2[wh2]))
if(l2>0) psub <- paste(sp2[wh2],sp2[wh2+1],sep=",") else psub <- sp2[wh2]
wh3 <- grep("))",psub)
el3 <- length(wh3)
if(el3>0) sp2[wh2] <- paste("param=c(",initrandomprec[1],",",initrandomprec[2],"))",sep="") else 
sp2[wh2] <- paste("param=c(",initrandomprec[1],",",initrandomprec[2],")",sep="")
if(l2 >0) sp2 <- sp2[-(wh2+1)]
} else {
lsp2 <- length(sp2)
llsp2 <- nchar(sp2[lsp2])
sp2[lsp2] <- substr(sp2[lsp2],1,(llsp2-1))
sp2 <- c(sp2,paste("param=c(",initrandomprec[1],",",initrandomprec[2],"))",sep=""))
}
newf <- sp2[1]; for(i in 2:length(sp2))  newf <- paste(newf,sp2[i],sep=",")
sp[wh] <- newf
fnew <- sp[1]; if(length(sp)>1) {for(i in 2:length(sp)) fnew <- paste(fnew,sp[i],sep="+")}
formret <- formula(paste("y ~ ",fnew,sep=""))
return(formret)
} else {
    return(form)
}
}
