#!/usr/bin/Rscript
cat(format(Sys.time(), "%b %d %X"),'!!!!!!!!!!!!!!!!\tstarted\t!!!!!!!!!!!!!!!!\n\n')
library(dplyr)
library(biomformat)
library(RJSONIO)
library(stringr)
library(data.table)
library(DBI)
library(MonetDB.R)
library(asarDB)
cat(format(Sys.time(), "%b %d %X"),'libraries loaded\n')
abundStatQ<-paste(" insert into asarLite.mv_abund_stat(sampleid,cnt,abundance)",
                  " select sampleid,",
                  " count(distinct m.md5id) as cnt,",
                  " sum(abundance) as abundance ",
                  " from abundance m",
                  " where sampleid in (select id as sampleid from smpl ",
                  "  except ",
                  "  select sampleid from mv_abund_stat)",
                  " group by sampleid")
addMD5Q<-'INSERT INTO md5(md5sum) values(?)'
dbname = "asar"
usr='asar'
pwd='asar'
webkey='y8MLy5NRMBfCRGzPcxYP9b86c'
#dtPath = '/Volumes/AS_WD_HFS/OIST/DBData/' #'/var/workspaceR/scalpelData/data/'
#tmpPath = '/Volumes/AS_WD_HFS/OIST/DBData/' #'/Volumes/AS_WD_HFS/Projects/R/MFCMAP/bash/' #'/var/workspaceR/scalpelData/data/'
rmdPath = '~/Documents/Projects/R/asarDB/load2DB.Rmd'
conn <- dbConnect(MonetDB.R::MonetDB.R(),
                  dbname = dbname,
                  user=usr,password=pwd,
                  timeout = 1200)
cat(format(Sys.time(), "%b %d %X"),'DB connected\n')
md5DB<-dbReadTable(conn,'md5')
cat(format(Sys.time(), "%b %d %X"),'MD5 loaded',dim(md5DB),'\n')
smpl<-dbReadTable(conn,'smpl')
cat(format(Sys.time(), "%b %d %X"),'Samples loaded',dim(smpl),'\n')
studyDB<-dbReadTable(conn,'study')
cat(format(Sys.time(), "%b %d %X"),'Projects loaded',dim(studyDB),'\n')
wd<-getwd()
fl<-dir(wd,pattern='project.mgp[0-9]+$')
for(f in fl){
  cat(format(Sys.time(), "%b %d %X"),f,'started\n')
  fdir<-normalizePath(f)
  setwd(fdir)
  path<-fdir
  proj.ID<-sub('project.','',f)
  idx<-which(studyDB$mgrastid==as.character(proj.ID))
  if(length(idx)==0){
    cat(format(Sys.time(), "%b %d %X"),'Metagenome project',as.character(proj.ID),' not found in the database\n')
    if((!exists('prjName'))){prjName<-proj.ID}
    if(!exists('prjDesc')){prjDesc<-paste0('project.',as.character(proj.ID))}
    if(nchar(trimws(prjName))==0){prjName<-as.character(proj.ID)}
    idStudy<-insertStudy2DB(conn,as.character(proj.ID),prjName,prjDesc)
    studyDB<-dbReadTable(conn,'study')
  }
  if(!file.exists('./abund.tsv')){
    idStudy<-studyDB$id[idx]
    mdt<-read.csv(file = paste0(proj.ID,'.meta.csv'))

    for( i in 1:dim(mdt)[1]){
      mgrastid<- mdt$MG.RAST.ID[i]
      cat(format(Sys.time(), "%b %d %X"),i,'Processing of ',as.character(mgrastid),' started\n')
      tbls<-dbGetQuery(conn,'SELECT * FROM smpl WHERE mgrastid = ?',mgrastid)
      if(dim(tbls)[1]==0){
        insertSample2DB(conn,path,idStudy,mgrastid,mdt$Metagenome.Name[i])
        cat(format(Sys.time(), "%b %d %X"),i,'Record of ',as.character(mgrastid),'(',smplid,') inserted into DB\n')
        tbls<-dbGetQuery(conn,'SELECT * FROM smpl WHERE mgrastid = ?',mgrastid)
      }
      smplid<-tbls$id
      cat(format(Sys.time(), "%b %d %X"),i,'Read of ',as.character(mgrastid),'(',smplid,') abundances started\n')
      rseed.fname<-normalizePath(gsub('//','/',paste0(path,'/seed',mgrastid,'.tsv')))
      fseed.fname<-normalizePath(gsub('//','/',paste0(path,'/fsub',mgrastid,'.tsv')))
      kegg.fname<-normalizePath(gsub('//','/',paste0(path,'/kegg',mgrastid,'.tsv')))
      if(file.exists(rseed.fname)&file.exists(fseed.fname)&file.exists(kegg.fname)){
        cat(format(Sys.time(), "%b %d %X"),i,'All files of ',as.character(mgrastid),'(',smplid,') found\n')
        rseed<-fread(rseed.fname,sep='$')
        names(rseed)<-c('md5','ab','usp')
        rseed$smplid<-smplid
        md5idx<-match(rseed$md5,md5DB$md5sum)
        rseed.len<-length(which(is.na(md5idx)))
        if(rseed.len>0){
          missingMD5<-unique(rseed$md5[which(is.na(md5idx))])
          for(s in missingMD5){rs<-dbSendStatement(conn,addMD5Q,s)}
          md5DB<-dbReadTable(conn,'md5')
          md5idx<-match(rseed$md5,md5DB$md5sum)
          rseed.len<-length(which(is.na(md5idx)))
        }
        rseed$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]

        fseed<-fread(fseed.fname,sep='$')
        names(fseed)<-c('md5','ab','ufun')
        fseed$smplid<-smplid
        md5idx<-match(fseed$md5,md5DB$md5sum)
        fseed.len<-length(which(is.na(md5idx)))
        if(fseed.len>0){
          missingMD5<-unique(fseed$md5[which(is.na(md5idx))])
          for(s in missingMD5){rs<-dbSendStatement(conn,addMD5Q,s)}
          md5DB<-dbReadTable(conn,'md5')
          md5idx<-match(fseed$md5,md5DB$md5sum)
          fseed.len<-length(which(is.na(md5idx)))
        }
        fseed$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]

        rkegg<-fread(kegg.fname,sep='$')
        names(rkegg)<-c('mgrastid','md5','ab','KO')
        rkegg$smplid<-smplid
        md5idx<-match(rkegg$md5,md5DB$md5sum)
        rkegg.len<-length(which(is.na(md5idx)))
        if(rkegg.len>0){
          missingMD5<-unique(rkegg$md5[which(is.na(md5idx))])
          for(s in missingMD5){rs<-dbSendStatement(conn,addMD5Q,s)}
          md5DB<-dbReadTable(conn,'md5')
          md5idx<-match(rkegg$md5,md5DB$md5sum)
          rkegg.len<-length(which(is.na(md5idx)))
        }
        rkegg$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]

        if((rkegg.len+fseed.len+rseed.len)==0){
          abund<-unique(rbind(rseed[,.(smplid,ab,md5id)],fseed[,.(smplid,ab,md5id)],rkegg[,.(smplid,ab,md5id)]))
          ab<-abund[,.(N=length(md5id),s=sum(ab),nab=length(unique(ab))),by=.(smplid,md5id)]
          if(max(ab$nab)==1){
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances started\n')
            fwrite(abund,file='./abund.tsv',col.names=FALSE,append=TRUE,quote=FALSE,sep='$',row.names=FALSE,na='')
            rm(abund,ab,rseed,fseed,rkegg)
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances finishes\n')
          }else{
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances impossible',dim(ab[nab>1])[1],' duplicates \n')
          }
        }else{
          cat(format(Sys.time(), "%b %d %X"),'Missing MD5s in ',as.character(mgrastid),'SEED:',rseed.len,'FSUB:',fseed.len,'KEGG:',rkegg.len,' \n')
        }
        cat(format(Sys.time(), "%b %d %X"),'Read of ',as.character(mgrastid),'(',i,') abundances finished\n')
      }else{
        cat(format(Sys.time(), "%b %d %X"),'Missing TSV files in ',as.character(mgrastid),' \n')
      }
    }
    afname<-normalizePath('./abund.tsv')
    if(file.exists(afname)){
      cat(format(Sys.time(), "%b %d %X"),'load',as.character(proj.ID),'abundances started\n')
      system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.abundance  FROM ',"'",afname,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
      cat(format(Sys.time(), "%b %d %X"),'load',as.character(proj.ID),'abundances finishes\n')
      rs<-dbSendStatement(conn,abundStatQ)
      cat(format(Sys.time(), "%b %d %X"),'update of ',as.character(proj.ID),'abundances stat finishes\n')
    }
  }else{
    cat(format(Sys.time(), "%b %d %X"),'Metagenome project',as.character(proj.ID),' have "abund.tsv" file\n')
  }
  setwd(wd)
  cat(format(Sys.time(), "%b %d %X"),f,'finished\n')
}