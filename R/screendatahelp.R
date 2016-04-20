screendatahelp <-
function(datset,variable,np = TRUE,ncpus=2,...){
#datset <- methdat;variable<- "groupfac"; np = TRUE;updateby=10000;ncpus=11
vari <- get(variable)
refu <- function(tp,...){
#varia <- vari
ifac <- is.factor(vari)
formch <- "y ~ x"
form <- y~x
el <- nrow(datset)
if(!ifac) print("Screen function only applies to factor variables") else {
    if(ncpus==1) {
        nlev <- length(levels(vari))
        if(nlev==2){
             if(np) pv <- apply(datset,1,function(vec) wilcox.test(form,data = cbind(y=as.numeric(vec),x=vari),exact=FALSE)[[3]]) else 
                    pv <- apply(datset,1,function(vec) t.test(form,data = cbind(y=as.numeric(vec),x=vari))[[3]])
                    } else {
             if(np) pv <- apply(datset,1,function(vec) kruskal.test(form,data = cbind(y=as.numeric(vec),x=vari))[[3]]) else 
                    pv <- apply(datset,1,function(vec) anova(lm(form,data = data.frame(y=as.numeric(vec),x=vari)))[[5]][1])     
                    }
             
        } else {
        #sfInit(parallel=TRUE,cpus=ncpus) 
        print("Exporting to slaves")
        sfExport("datset","vari","formch")
        print("Exporting done, start computing p-values")
        nlev <- length(levels(vari))
        if(nlev==2){
             if(np) pv <- sfApply(datset,1,function(vec) wilcox.test(as.formula(formch),data = cbind(y=as.numeric(vec),x=vari),exact=FALSE)[[3]]) else 
                    pv <- sfApply(datset,1,function(vec) t.test(as.formula(formch),data = cbind(y=as.numeric(vec),x=vari))[[3]])
                    } else {
             if(np) pv <- sfApply(datset,1,function(vec) kruskal.test(as.formula(formch),data = cbind(y=as.numeric(vec),x=vari))[[3]]) else 
                    pv <- sfApply(datset,1,function(vec) anova(lm(as.formula(formch),data = data.frame(y=as.numeric(vec),x=vari)))[[5]][1])     
                    }
        sfRemoveAll()
        } 
        return(pv)
        }}
    pval <- refu(1)
    return(pval)
}

