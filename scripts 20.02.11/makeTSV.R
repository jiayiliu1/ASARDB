#!/usr/bin/Rscript
args <- commandArgs(trailingOnly = TRUE)
print(args)
proj.ID<-args[1]
library(data.table)
library(stringr)
library(plyr)
library(biomformat)
library(asarDB)
mdt<-read.csv(file = paste0(proj.ID,'.meta.csv'))
i<-1
path<-'./'
for(mgrastid in mdt$MG.RAST.ID){
  rseed.fname<-normalizePath(gsub('//','/',paste0(path,'/seed',mgrastid,'.tsv')))
  cat(format(Sys.time(), "%b %d %X"),'loadSeed',mgrastid,rseed.fname,'\n')
  if(!file.exists(rseed.fname)){
    s<-loadMGRAST(path,paste0(mgrastid,'.seed'))
    sa<-unique(s[,.(ab=length(`query sequence id`),`semicolon separated list of annotations`),by=.(`hit m5nr id (md5sum)`)])
    sa.1<-sa[,list(ab,sp=unlist(gsub('(\\]|\\[)','',unlist( strsplit(`semicolon separated list of annotations`, "\\]; *\\[" ) )))),
             by=.(`hit m5nr id (md5sum)`)]
    cat(format(Sys.time(), "%b %d %X"),'readSEED',mgrastid,'\n')
    fwrite(sa.1,file=rseed.fname,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
  }
  fseed.fname<-normalizePath(gsub('//','/',paste0(path,'/fsub',mgrastid,'.tsv')))
  cat(format(Sys.time(), "%b %d %X"),'loadFSeed',mgrastid,fseed.fname,'\n')
  if(!file.exists(fseed.fname)){
    f<-loadMGRAST(path,paste0(mgrastid,'.fsub'))
    fa<-unique(f[,.(ab=length(`query sequence id`),`semicolon separated list of annotations`),by=.(`hit m5nr id (md5sum)`)])
    fa.1<-fa[,list(ab,fun=unlist(gsub('accession=\\[SS([0-9]+)\\].*','SS\\1',
                                      unlist(str_split(`semicolon separated list of annotations`,';'))))),
             by=.(`hit m5nr id (md5sum)`)]
    cat(format(Sys.time(), "%b %d %X"),'readSEED',mgrastid,'\n')
    fwrite(fa.1,file=fseed.fname,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
  }
  kegg.fname<-normalizePath(gsub('//','/',paste0(path,'/kegg',mgrastid,'.tsv')))
  cat(format(Sys.time(), "%b %d %X"),'loadKEGG',mgrastid,kegg.fname,'\n')
  if(!file.exists(kegg.fname)){
    rkegg<-readKEGG(mgrastid,path)
    cat(format(Sys.time(), "%b %d %X"),'readKEGG',mgrastid,'\n')
    writeDTannot(rannot = rkegg,mgrastID = mgrastid,fname = kegg.fname,annot.col = 'ko')
  }
}
