nms<-read.delim('~/Downloads/new_taxdump/names.dmp',header=FALSE,sep='|')
names(nms)<-c('taxid','name','uname','nclass')
nms$name<-gsub('\t','',nms$name)
nms$uname<-gsub('\t','',nms$uname)
nms$nclass<-gsub('\t','',nms$nclass)
nms<-nms[,-5]

fl<-dir(path='~/Documents/',pattern='pathview.*.Rdata',recursive=TRUE)
cat(format(Sys.time(), "%b %d %X"),dim(taxDF),dim(funDF),'\n')
taxL<-list()
funL<-list()
for(f in fl){
load(paste0('~/Documents/',f))
taxL[[length(taxL)+1]]<-unique(funtaxall[,c("usp","species","genus","family","order","class","phylum","domain")])
funL[[length(funL)+1]]<-unique(funtaxall[,c("ufun","FUN2","FUN3","FUN4")])
cat(format(Sys.time(), "%b %d %X"),f,length(taxL),length(funL),'\n')
}
taxDF<-unique(do.call(rbind,taxL))
funDF<-unique(do.call(rbind,funL))
cat(format(Sys.time(), "%b %d %X"),f,dim(taxDF),dim(funDF),'\n')
save(fl,taxDF,funDF,f,file='dbprep.RData')
flW<-dir(path='/Volumes/AS_WD_HFS/',pattern='pathview.*.Rdata',recursive=TRUE)
cat(format(Sys.time(), "%b %d %X"),dim(taxDF),dim(funDF),'\n')
for(f in flW){
load(paste0('/Volumes/AS_WD_HFS/',f))
taxL[[length(taxL)+1]]<-unique(funtaxall[,c("usp","species","genus","family","order","class","phylum","domain")])
funL[[length(funL)+1]]<-unique(funtaxall[,c("ufun","FUN2","FUN3","FUN4")])
cat(format(Sys.time(), "%b %d %X"),f,length(taxL),length(funL),'\n')
}
taxDF<-unique(do.call(rbind,taxL))
funDF<-unique(do.call(rbind,funL))
save(fl,taxDF,funDF,f,nms,file='dbprep.RData')
usps<-unique(as.character(taxDF$usp))
uspDF<-data.frame(usp=usps,name=usps,taxid=-1)
uspDF$name<-gsub('^(Siphoviridae|Myoviridae|Caudovirales|44AHJD-like phages|Viruses|Tectivirus|Autographivirinae|Inovirus|Podoviridae) *','',uspDF$name)
uspDF$name<-gsub("^.*viruses *","",uspDF$name)
uspDF$name<-gsub("(viruses|Viruses) *","",uspDF$name)
uspDF$name<-gsub('plasmid.*$','',uspDF$name)
uspDF$name<-gsub("'+","'",uspDF$name)
uspDF$name<-gsub(" +\\.$","",uspDF$name)
uspDF$name<-trimws(uspDF$name)
idx<-match(uspDF$usp,nms$name)
uspDF$taxid[which(!is.na(idx))]<-nms$taxid[idx[which(!is.na(idx))]]
idxNF<-which(uspDF$taxid<0)
idx<-match(uspDF$name[idxNF],nms$name)
uspDF$taxid[idxNF[which(!is.na(idx))]]<-nms$taxid[idx[which(!is.na(idx))]]
idxNF<-which(uspDF$taxid<0)
uspDF$name2<-uspDF$name
uspDF$name2[idxNF]<-gsub(' +[^ ]+$','',uspDF$name[idxNF])
idx<-match(uspDF$name2[idxNF],nms$name)
uspDF$taxid[idxNF[which(!is.na(idx))]]<-nms$taxid[idx[which(!is.na(idx))]]
idxNF<-which(uspDF$taxid<0)
write.table(uspDF[idxNF,],file='unfound4.txt',sep='\t')
save(uspDF,file='taxInit.RData')
uspDFf<-uspDF[uspDF$taxid>0,]
 resSpCLS<-list()
 splSp<-split(uspDFf,c(rep(1:(dim(uspDFf)[1]/50),50),1:(dim(uspDFf)[1]%%50)))
 i1<-length(resSpCLS)+1;i2<-length(splSp)
for(i in i1:i2){
 sub<-splSp[[i]]
 taxCLS<-classification(sub$taxid,db='ncbi')
 resSpCLS[[i]]<-taxCLS
 cat('sp',format(Sys.time(), "%b %d %X"),i,i2,'\n')
 if(i%%50==0){
 save(res,resSpCLS,splSp,uspDFf,file='taxRes.RData')
 cat('new save\n')
 }
 }
 save(res,resSpCLS,splSp,uspDFf,uspDF,file='taxInit.RData')
 resU<-unlist(resSpCLS,recursive=FALSE)
 resSp<-sapply(resU,function(.x)ifelse(any(.x$rank=='species'),.x$name[.x$rank=='species'],NA))
 resGen<-sapply(resU,function(.x)ifelse(any(.x$rank=='genus'),.x$name[.x$rank=='genus'],NA))
 resFun<-sapply(resU,function(.x)ifelse(any(.x$rank=='family'),.x$name[.x$rank=='family'],NA))
 resOrd<-sapply(resU,function(.x)ifelse(any(.x$rank=='order'),.x$name[.x$rank=='order'],NA))
 resCl<-sapply(resU,function(.x)ifelse(any(.x$rank=='class'),.x$name[.x$rank=='class'],NA))
 resPh<-sapply(resU,function(.x)ifelse(any(.x$rank=='phylum'),.x$name[.x$rank=='phylum'],NA))
 resDom<-sapply(resU,function(.x)ifelse(any(.x$rank=='superkingdom'),.x$name[.x$rank=='superkingdom'],NA))
 taxDF<-uspDFf[,c('taxid','usp','name')]
 taxDF$species<-resSp[match(taxDF$taxid,names(resSp))]
 taxDF$genus<-resGen[match(taxDF$taxid,names(resGen))]
 taxDF$family<-resFun[match(taxDF$taxid,names(resFun))]
 taxDF$ordr<-resOrd[match(taxDF$taxid,names(resOrd))]
 taxDF$class<-resCl[match(taxDF$taxid,names(resCl))]
 taxDF$phylum<-resPh[match(taxDF$taxid,names(resPh))]
 taxDF$dom<-resDom[match(taxDF$taxid,names(resDom))]
save(res,resSpCLS,splSp,uspDFf,uspDF,taxDF,file='taxInit.RData')