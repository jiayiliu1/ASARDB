cat(format(Sys.time(), "%b %d %X"),'!!!!!!!!!!!!!!!!\tstarted\t!!!!!!!!!!!!!!!!\n\n')
library(dplyr)
library(biomformat)
library(RJSONIO)
library(stringr)
library(data.table)
library(DBI)
library(MonetDBLite)
library(asarDB)
cat(format(Sys.time(), "%b %d %X"),'libraries loaded\n')
dbname = "asar"
usr='asar'
pwd='asar'
webkey='y8MLy5NRMBfCRGzPcxYP9b86c'
dtPath = '/Volumes/AS_WD_HFS/OIST/DBData/' #'/var/workspaceR/scalpelData/data/'
tmpPath = '/Volumes/AS_WD_HFS/OIST/DBData/' #'/Volumes/AS_WD_HFS/Projects/R/MFCMAP/bash/' #'/var/workspaceR/scalpelData/data/'
rmdPath = '~/Documents/Projects/R/asarDB/load2DB.Rmd'
conn <- dbConnect(MonetDBLite::MonetDBR(),
                  dbname = dbname,
                  user=usr,password=pwd,
                  timeout = 1200)
makeMD5<-function(rseedN,fseedN,keggN){

}
cat(format(Sys.time(), "%b %d %X"),'DB connected\n')
md5DB<-dbReadTable(conn,'md5')
cat(format(Sys.time(), "%b %d %X"),'MD5 loaded\n')
smpl<-dbReadTable(conn,'smpl')
cat(format(Sys.time(), "%b %d %X"),'Samples loaded\n')
studyDB<-dbReadTable(conn,'study')
cat(format(Sys.time(), "%b %d %X"),'Projects loaded\n')
####### KOrts #######
# list of KOs already in the database
kort<-dbReadTable(conn,'kort')
cat(format(Sys.time(), "%b %d %X"),'KEGG Orthology loaded\n')
####### DB samples #######
# list of samples already in the database
loadedSmpl<-dbGetQuery(conn,'select mgrastid from smpl s join (select distinct sampleid from abundance) i on i.sampleid=s.id')
cat(format(Sys.time(), "%b %d %X"),'Processed samples loaded\n')
####### SEED #######
# list of SEED IDs already in the database
seedLite<-dbGetQuery(conn,'select id,accession from seedlite')
cat(format(Sys.time(), "%b %d %X"),'SEED loaded\n')
####### Taxonomy #######
# list of taxonomy already in the database.
# Match between MG-RAST and taxonomy by matching taxonomy name and MG-RAST usp.
# At the moment in the taxonomy
taxLite<-dbGetQuery(conn,'select id,usp,name from taxlite')
cat(format(Sys.time(), "%b %d %X"),'Taxonomy is loaded\n')
wd<-getwd()
missingAbundFile<-paste0(wd,'/missing.abund.tsv')
missingMD5File<-paste0(wd,'/missing.md5.tsv')
missingKEGGFile<-paste0(wd,'/missing.KEGG.tsv')
missingSLFile<-paste0(wd,'/missing.SL.tsv')
missingTaxFile<-paste0(wd,'/missing.TAX.tsv')
fl<-dir(wd,pattern='project.mgp[0-9]+$')
missingAbund<-list()
missingKO<-list()
missingSL<-list()
missingTAX<-list()
######### Main cycle over project folders #######
for(f in fl){
  cat(format(Sys.time(), "%b %d %X"),f,'started\n')
  fdir<-normalizePath(paste0(wd,'/',f))
  setwd(fdir)
  path<-fdir
  proj.ID<-sub('project.','',f)
  idx<-which(studyDB$mgrastid==proj.ID)
  if(length(idx)>0){
    ########### Project exists in the DB ##########
    idStudy<-studyDB$id[idx]
    mdt<-read.csv(file = paste0(proj.ID,'.meta.csv'))
    idx<-match(as.character(mdt$MG.RAST.ID),loadedSmpl$mgrastid)
    if(length(which(is.na(idx)))>0){
      ############# There are unloaded metagenomes in the project ############
      mgids<-which(is.na(idx))
      for(i in mgids){
        mgrastid<-as.character(mdt$MG.RAST.ID[i])
        tbls<-dbGetQuery(conn,'SELECT * FROM smpl WHERE mgrastid = ?',mgrastid)
        if(dim(tbls)[1]==0){
          insertSample2DB(conn,path,idStudy,mgrastid,mdt$Metagenome.Name[i])
          tbls<-dbGetQuery(conn,'SELECT * FROM smpl WHERE mgrastid = ?',mgrastid)
        }
        smplid<-tbls$id
        cat(format(Sys.time(), "%b %d %X"),i,'Read of ',as.character(mgrastid),'(',smplid,') abundances started\n')
        rseed.fname<-normalizePath(gsub('//','/',paste0(path,'/seed',mgrastid,'.tsv')))
        fseed.fname<-normalizePath(gsub('//','/',paste0(path,'/fsub',mgrastid,'.tsv')))
        kegg.fname<-normalizePath(gsub('//','/',paste0(path,'/kegg',mgrastid,'.tsv')))
        if(file.exists(rseed.fname)&file.exists(fseed.fname)&file.exists(kegg.fname)){
          ########### All files are in the folder #########
          rseed<-fread(rseed.fname,sep='$')
          names(rseed)<-c('md5','ab','usp')
          rseed$smplid<-smplid
          md5idx<-match(rseed$md5,md5DB$md5sum)
          rseed.len<-length(which(is.na(md5idx)))
#          rseed$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
          if(rseed.len>0){
            ########## Save missing Taxonomy md5 for further processing #########
            newTax<-rseed[is.na(md5idx)]
            missingTAX[[length(missingTAX)+1]]<-newTax
          }

          fseed<-fread(fseed.fname,sep='$')
          names(fseed)<-c('md5','ab','ufun')
          fseed$smplid<-smplid
          md5idx<-match(fseed$md5,md5DB$md5sum)
          fseed.len<-length(which(is.na(md5idx)))
#          fseed$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
          if(fseed.len>0){
            ########## Save missing SEED md5 for further processing #########
            newSL<-fseed[is.na(md5idx)]
            missingSL[[length(missingSL)+1]]<-newSL
          }

          rkegg<-fread(kegg.fname,sep='$')
          names(rkegg)<-c('mgrastid','md5','ab','KO')
          rkegg$smplid<-smplid
          md5idx<-match(rkegg$md5,md5DB$md5sum)
          rkegg.len<-length(which(is.na(md5idx)))
#          rkegg$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
          md5s<-c()
          if((rkegg.len)>0){
            ########## Save missing KO md5 for further processing #########
           newKegg<-rkegg[is.na(md5idx)]
            missingKO[[length(missingKO)+1]]<-newKegg
          }
          abund<-unique(rbind(rseed[,.(smplid,ab,md5)],fseed[,.(smplid,ab,md5)],rkegg[,.(smplid,ab,md5)]))
          ab<-abund[,.(N=length(md5),s=sum(ab),nab=length(unique(ab))),by=.(smplid,md5)]
          if(max(ab$nab)==1){
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances started\n')
            missingAbund[[length(missingAbund)+1]]<-abund
            rm(abund,ab,rseed,fseed,rkegg)
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances finishes\n')
          }else{
            cat(format(Sys.time(), "%b %d %X"),'write',as.character(mgrastid),'abundances impossible',dim(ab[nab>1])[1],' duplicates \n')
          }
        }else{
          cat(format(Sys.time(), "%b %d %X"),'Missing TSV files in ',as.character(mgrastid),' \n')
        }
      }
    }else{
      cat(format(Sys.time(), "%b %d %X"),'Metagenome project',as.character(proj.ID),' completely loaded to the database\n')
    }
  }else{
    cat(format(Sys.time(), "%b %d %X"),'Metagenome project',as.character(proj.ID),' not found in the database\n')
  }
  setwd(wd)
  cat(format(Sys.time(), "%b %d %X"),f,'finished\n')
}

mKO<-do.call(rbind,missingKO)
newKO<-unique(mKO[,.(md5,KO)])
mSL<-do.call(rbind,missingSL)
newSL<-unique(mSL[,.(md5,ufun)])
mTax<-do.call(rbind,missingTAX)
newTax<-unique(mTax[,.(md5,usp)])
newMD5<-unique(rbind(newTax[,.(md5)],newSL[,.(md5)],newKO[,.(md5)]))
cat(format(Sys.time(), "%b %d %X"),' Aggregate abundances starts\n')
mAbund<-do.call(rbind,missingAbund)
ab<-mAbund[,.(N=length(md5),s=sum(ab),nab=length(unique(ab))),by=.(smplid,md5)]
cat(format(Sys.time(), "%b %d %X"),' Aggregate abundances done\n')
idx<-match(newMD5$md5,md5DB$md5sum)
if(any(!is.na(idx))){
  cat(format(Sys.time(), "%b %d %X"),' Not all new MD5 is missing from DB\n')
}else{
  ######## Processing missing MD5 ########
  cat(format(Sys.time(), "%b %d %X"),' Update MD5s starts.\n')
  fwrite(newMD5,file=missingMD5File,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
  system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.md5(md5sum)  FROM ',"'",missingMD5File,"'(md5sum)",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
  cat(format(Sys.time(), "%b %d %X"),' Update MD5s done.\n')
  md5DB<-dbReadTable(conn,'md5')
  cat(format(Sys.time(), "%b %d %X"),'MD5s reloaded from DB.\n')
  md5idx<-match(mAbund$md5,md5DB$md5sum)
  mAbund.len<-length(which(is.na(md5idx)))
  if(mAbund.len==0){
    cat(format(Sys.time(), "%b %d %X"),' Update abundances starts.\n')
    mAbund$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
    fwrite(mAbund[,.(smplid,ab,md5id)],file=missingAbundFile,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
    system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.abundance  FROM ',"'",missingAbundFile,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
    cat(format(Sys.time(), "%b %d %X"),' Update abundances done.\n')
  }else{
    cat(format(Sys.time(), "%b %d %X"),mAbund.len,' MD5s is still missing from DB.\n')
  }
  ######## Taxonomy reload ########
  md5idx<-match(newTax$md5,md5DB$md5sum)
  newTax$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
  newTax.len<-length(which(is.na(md5idx)))
  if(newTax.len==0){
    cat(format(Sys.time(), "%b %d %X"),' Update taxonomy starts.\n')
    taxidx<-match(newTax$usp,taxLite$name)
    if(any(is.na(taxidx))){
      newTax.names<-unique(newTax$usp[is.na(taxidx)])
      cat(format(Sys.time(), "%b %d %X"),' Load new ',length(newTax),' taxonomy starts.\n')
      for(k in newTax.names){
        rs<-dbSendStatement(conn,
                            paste('INSERT INTO',
                                  'taxlite (usp,name)',
                                  'VALUES (?,?);'),
                            k,k)
        dbHasCompleted(rs)
        dbGetRowsAffected(rs)
        dbClearResult(rs)
        cat(format(Sys.time(), "%b %d %X"),k,' Update taxonomy.\n')
      }
      taxLite<-dbGetQuery(conn,'select id,usp,name from taxlite')
      cat(format(Sys.time(), "%b %d %X"),'Taxonomy is reloaded\n')
      taxidx<-match(newTax$usp,taxLite$name)
    }
    newTax$taxid<-taxLite$id[taxidx]
    fwrite(unique(newTax[,.(taxid,md5id)]),file=missingTaxFile,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
    system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.md5_taxlite  FROM ',"'",missingTaxFile,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
    cat(format(Sys.time(), "%b %d %X"),' Update taxonomy done.\n')
  }else{
    cat(format(Sys.time(), "%b %d %X"),newTax.len,' MD5s is still missing from DB.\n')
  }
  ######## SEED reload ########
  md5idx<-match(newSL$md5,md5DB$md5sum)
  newSL$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
  slidx<-match(newSL$ufun,seedLite$accession)
   if(any(is.na(slidx))){
     newSL.ufun<-uniqu(newSL$ufun[is.na(slidx)])
     for(k in newSL.ufun){
      rs<-dbSendStatement(conn,
                          paste('INSERT INTO',
                                'seedlite (accession,fun1,fun2,fun3,fun4)',
                                'VALUES (?,?,?,?,?);'),
                          k,'-1','-1','-1','-1')
      dbHasCompleted(rs)
      dbGetRowsAffected(rs)
      dbClearResult(rs)
      cat(format(Sys.time(), "%b %d %X"),k,' Update SEED \n')
     }
     seedLite<-dbGetQuery(conn,'select id,accession from seedlite')
     cat(format(Sys.time(), "%b %d %X"),'SEED loaded\n')
     slidx<-match(newSL$ufun,seedLite$accession)
   }
   #careful examination is needed as all nonmatching ufuns happened to be strings like
   #Outer membrane receptor for ferrienterochelin and colicins]
  newSL$slid<-seedLite$id[slidx]
  fwrite(unique(newSL[,.(slid,md5id)]),file=missingSLFile,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
  system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.seedlite_md5  FROM ',"'",missingSLFile,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
  cat(format(Sys.time(), "%b %d %X"),' Update SEED done.\n')

  md5idx<-match(newKO$md5,md5DB$md5sum)
  newKO$md5id[!is.na(md5idx)]<-md5DB$id[md5idx[!is.na(md5idx)]]
  koidx<-match(newKO$KO,kort$access)
  if(any(is.na(koidx))){
     newKO.access<-newKO$KO[is.na(koidx)]
    for(k in newKO.access){
      rs<-dbSendStatement(conn,
                          paste('INSERT INTO',
                                'kort (access,funct,description)',
                                'VALUES (?,?,?);'),
                          k,'-1','-1')
      dbHasCompleted(rs)
      dbGetRowsAffected(rs)
      dbClearResult(rs)
    }
    kort<-dbReadTable(conn,'kort')
    cat(format(Sys.time(), "%b %d %X"),'KEGG Orthology loaded\n')
    koidx<-match(newKO$KO,kort$access)
  }
  newKO$koid<-kort$id[koidx]
  fwrite(unique(newKO[,.(koid,md5id)]),file=missingKEGGFile,col.names=FALSE,append=FALSE,quote=FALSE,sep='$',row.names=FALSE,na='')
  system.time(system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.md5_kort  FROM ',"'",missingKEGGFile,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut)
  cat(format(Sys.time(), "%b %d %X"),' Update KEGG done.\n')

}
# md5idx<-match(rseed$md5,md5DB$md5sum)
# md5s<-c(md5s,rseed$md5[is.na(md5idx)])
# taxidx<-match(rseed$usp,taxLite$usp)
# if(any(is.na(slidx))){
#   newTax<-rseed$ufun[is.na(slidx)]
#   for(k in newSL){
#     rs<-dbSendStatement(conn,
#                         paste('INSERT INTO',
#                               'taxlite (usp,name)',
#                               'VALUES (?,?,?,?);'),
#                         k,k)
#     dbHasCompleted(rs)
#     dbGetRowsAffected(rs)
#     dbClearResult(rs)
#   }
#   seedLite<-dbGetQuery(conn,'select id,accession from seedlite')
#   cat(format(Sys.time(), "%b %d %X"),'SEED loaded\n')
#   slidx<-match(rseed$ufun,seedLite$accession)
# }
# rseed$slid<-seedLite$id[slidx]


# md5idx<-match(fseed$md5,md5DB$md5sum)
# md5s<-c(md5s,fseed$md5[is.na(md5idx)])
# slidx<-match(fseed$ufun,seedLite$accession)
# if(any(is.na(slidx))){
#   newSL<-fseed$ufun[is.na(slidx)]
#   for(k in newSL){
#     rs<-dbSendStatement(conn,
#                         paste('INSERT INTO',
#                               'seedlite (accession,fun1,fun2,fun3,fun4)',
#                               'VALUES (?,?,?,?);'),
#                         k,'prjID','prjDesc','-1')
#     dbHasCompleted(rs)
#     dbGetRowsAffected(rs)
#     dbClearResult(rs)
#   }
#   seedLite<-dbGetQuery(conn,'select id,accession from seedlite')
#   cat(format(Sys.time(), "%b %d %X"),'SEED loaded\n')
#   slidx<-match(fseed$ufun,seedLite$accession)
# }
# fseed$slid<-seedLite$id[slidx]

# md5idx<-match(rkegg$md5,md5DB$md5sum)
# md5s<-c(md5s,rkegg$md5[is.na(md5idx)])
# koidx<-match(rkegg$KO,kort$access)
# if(any(is.na(koidx))){
#   newKO<-rkegg$KO[is.na(koidx)]
#   for(k in newKO){
#     rs<-dbSendStatement(conn,
#                         paste('INSERT INTO',
#                               'kort (access,funct,description)',
#                               'VALUES (?,?,?);'),
#                         k,'prjID','prjDesc')
#     dbHasCompleted(rs)
#     dbGetRowsAffected(rs)
#     dbClearResult(rs)
#   }
#   kort<-dbReadTable(conn,'kort')
#   cat(format(Sys.time(), "%b %d %X"),'KEGG Orthology loaded\n')
#   koidx<-match(rkegg$KO,kort$access)
# }
# rkegg$koid<-kort$id[koidx]
