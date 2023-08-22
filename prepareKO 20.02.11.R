library(KEGGREST)
load('~/Documents/Projects/R/asarDB/data/keggmappings.Rdata')
kortIDs<-unique(kegg$K)
spl<-split(kortIDs,c(rep(1:(length(kortIDs)/9),9),1:(length(kortIDs)%%9)))
resKO<-list()
i1<-length(resKO)+1;i2<-length(spl)
for(i in i1:i2){
 set<-spl[[i]]
 resKO[[i]]<-keggGet(set)
 cat('ko',format(Sys.time(), "%b %d %X"),i,i2,'\n')
 if(i%%50==0){
   save(resKO,spl,kortIDs,file='koRes.RData')
   cat('new save\n')
 }
}
save(resKO,spl,kortIDS,file='koRes.RData')

resU<-unlist(resKO,recursive=FALSE)
access<-sapply(resU,function(.x).x$ENTRY)
name<-sapply(resU,function(.x).x$NAME)
descr<-sapply(resU,function(.x).x$DEFINITION)
kortDF<-data.frame(access=access,stringsAsFactors=FALSE)
kortDF$name<-name
kortDF$description<-descr
kortDF$id<-as.integer(sub('^K','',kortDF$access))
save(kortDF,resKO,spl,kortIDs,file='koRes.RData')