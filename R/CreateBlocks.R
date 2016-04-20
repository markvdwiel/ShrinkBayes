CreateBlocks <-
function(ds,size=10000){
nfeat <- nrow(ds)
sec <- seq(1,nfeat,by=size)
secl <- length(sec)
if(sec[secl] > nfeat-2) sec[secl] <- nfeat+1 else sec <- c(sec,nfeat+1)
blocks <- c()
for(i in 2:length(sec)) blocks <- rbind(blocks,c(sec[i-1],sec[i]-1))
return(blocks)
}

