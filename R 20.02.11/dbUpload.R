#' Function to write data.table with SEED annotation, i.e. mappings between md5,
#' function and taxonomy with abundance information in it.
#'
#' @param rseed data.table to get data from
#' @param mgrastID ID of the metagenome to write
#' @param fname where to write into
#'
#' @export
#'
writeDTseed<-function(rseed,mgrastID,fname){
  #cat('Writing scan',md$number,md$retentionTime,'\n')
  resDF<-data.table(mgrastid=mgrastID,
                    md5sum=rseed$md5sum,
                    ab=rseed$ab,
                    ufun=gsub('\\$',':',rseed$ufun),
                    usp=gsub('\\$',':',rseed$usp))
  resDF<-unique(resDF)
  if(file.exists(fname)){
    app=TRUE
  }else{
    app=FALSE
  }
  fwrite(resDF,file=fname,col.names=FALSE,append=app,quote=FALSE,sep='$',row.names=FALSE,na='')
}

#' Function to write annotation table, i.e. mappings between md5 and external db reference, which do not have
#' abundance information in it.
#'
#' @param rannot data.table with mapping. It should have column 'md5sum'
#' @param mgrastID string with mgrast ID of the metagenome
#' @param fname file to write data.table
#' @param annot.col annotation column name
#'
#' @export
#'
writeDTannot<-function(rannot,mgrastID,fname,annot.col='ko'){
  #cat('Writing scan',md$number,md$retentionTime,'\n')
  resDF<-data.table(mgrastid=mgrastID,
                    md5sum=rannot$md5sum,
                    ab=rannot$ab,
                    annot=gsub('\\$',':',as.data.frame(rannot[,annot.col, with=FALSE])[,1]))
  resDF<-unique(resDF)
  if(file.exists(fname)){
    app=TRUE
  }else{
    app=FALSE
  }
  fwrite(resDF,file=fname,col.names=FALSE,append=app,quote=FALSE,sep='$',row.names=FALSE,na='')
}

#' Using bulk load to populate the database
#'
#' @param fname file to load
#'
#' @export
#'
uploadSeedDT<-function(fname){
  sigFname<-normalizePath(fname)
  system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.tmp_seed  FROM ',"'",sigFname,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut
  cat(format(Sys.time(), "%b %d %X"),'\n',mcOut,'\n')
  system2('wc', sigFname,stdout=TRUE)->wc
  wcn<-as.numeric(unlist(strsplit(trimws(wc),' '))[1:3])
  cat(wcn[1],'\n')

}


uploadAnnotDT<-function(fname){
  sigFname<-normalizePath(fname)
  system(paste0('mclient  -d asar  -s  "COPY INTO asarlite.tmp_annot  FROM ',"'",sigFname,"'",' DELIMITERS \'$\';"'),intern = TRUE)->mcOut
  cat(format(Sys.time(), "%b %d %X"),'\n',mcOut,'\n')
  system2('wc', sigFname,stdout=TRUE)->wc
  wcn<-as.numeric(unlist(strsplit(trimws(wc),' '))[1:3])
  cat(wcn[1],'\n')

}

#' Load  SEED annotation for particular metagenome to the database
#'
#' @param conn connection to the DB
#' @param mgrastid metagenome id
#' @param path where intermediate files will be stored
#'
#' @export
#'
loadSeed2DB<-function(conn,mgrastid,path){
  rseed.fname<-normalizePath(gsub('//','/',paste0(path,'/rseed',mgrastid,'.tsv')))
  cat(format(Sys.time(), "%b %d %X"),'loadSeed',mgrastid,rseed.fname,'\n')
  if(!file.exists(rseed.fname)){
    rseed<-readSEED(mgrastid,path)
    cat(format(Sys.time(), "%b %d %X"),'readSEED',mgrastid,'\n')
    writeDTseed(rseed,mgrastid,rseed.fname)
  }
  beforeTMP<-dbGetQuery(conn,'select count(*) as num from tmp_seed;')
  if(beforeTMP$num>0){
    rs<-dbSendStatement(conn,'delete from tmp_seed;')
  }
  uploadSeedDT(rseed.fname)
  smpl<-dbGetQuery(conn,'select * from smpl where mgrastid=?',mgrastid)
  first<-dbGetQuery(conn,paste('SELECT count(distinct md5sum) as cnt FROM (',
                               '(SELECT md5sum FROM tmp_seed) EXCEPT (SELECT md5sum FROM asarLite.md5)',
                               ') not_in_md5'))

  newMD5<-first$cnt
  if(newMD5>0){
    newMD5s<-dbGetQuery(conn,paste('SELECT distinct tmp_seed.md5sum FROM (',
                                   '(SELECT md5sum FROM tmp_seed) EXCEPT (SELECT md5sum  FROM asarLite.md5)',
                                   ') not_in_md5 JOIN tmp_seed ON not_in_md5.md5sum=tmp_seed.md5sum'))
    rs<-dbSendStatement(conn,paste('insert into md5 (md5sum) (select distinct md5sum from (',
                                   '(SELECT distinct md5sum FROM tmp_seed) EXCEPT (SELECT md5sum  FROM asarLite.md5)',
                                   ') not_in_md5)'))
  }
  #dbWriteTable(conn,'mdtax_tmp',unique(rseed[,c('md5sum','usp')]),overwrite=TRUE,temporary=TRUE)
  newUSP<-dbGetQuery(conn,paste('SELECT distinct tmp_seed.usp FROM (',
                                '(SELECT usp FROM tmp_seed) EXCEPT (SELECT name as usp  FROM asarLite.taxlite)',
                                ') not_in_tax JOIN tmp_seed ON not_in_tax.usp=tmp_seed.usp'))
  if(dim(newUSP)[1]>0){
    cat('New USP found:',dim(newUSP)[1],'\n')
    uspDF<-data.frame(usp=newUSP$usp,name=newUSP$usp,taxid=-1)
    # uspDF$name<-gsub('^(Siphoviridae|Myoviridae|Caudovirales|44AHJD-like phages|Viruses|Tectivirus|Autographivirinae|Inovirus|Podoviridae) *','',uspDF$name)
    # uspDF$name<-gsub("^.*viruses *","",uspDF$name)
    # uspDF$name<-gsub("(viruses|Viruses) *","",uspDF$name)
    # uspDF$name<-gsub('plasmid.*$','',uspDF$name)
    uspDF$usp<-gsub("'+","'",uspDF$usp)
    uspDF$usp<-gsub(" +\\.$","",uspDF$usp)
    uspDF$usp<-trimws(uspDF$usp)
    dbWriteTable(conn,'newtax_tmp',uspDF,overwrite=TRUE,temporary=TRUE)
    rs<-dbSendStatement(conn,'insert into taxlite (taxid,usp,name) (select taxid,usp,name from newtax_tmp);')
  }
  newUSPmap<-dbGetQuery(conn,paste('select count(*) as cnt from (SELECT distinct tmp_seed.md5sum,tmp_seed.usp FROM (',
                                   '(SELECT md5sum,usp FROM tmp_seed) EXCEPT (SELECT d.md5sum as md5sum,t.name as usp FROM md5_taxLite m ',
                                   'join md5 d on m.md5id=d.id',
                                   'join taxlite t on m.taxLiteID=t.id)',
                                   ') not_in_map JOIN tmp_seed ON not_in_map.usp=tmp_seed.usp and not_in_map.md5sum=tmp_seed.md5sum) mdtax_cnt '))
  if(newUSPmap$cnt>0){
    rs<-dbSendStatement(conn,paste('insert into md5_taxlite (md5id,taxLiteID) ',
                                   '(SELECT distinct  p.md5id,x.id as taxLiteID from',
                                   '(SELECT distinct tmp_seed.md5sum,tmp_seed.usp FROM (',
                                   '(SELECT d.id as md5id,md5sum,usp FROM tmp_seed join md5 d on  tmp_seed.md5sum=d.md5sum) EXCEPT (SELECT md5id,d.md5sum as md5sum,t.name as usp FROM md5_taxLite m ',
                                   'join md5 d on m.md5id=d.id',
                                   'join taxlite t on m.taxLiteID=t.id)',
                                   ') not_in_map JOIN tmp_seed ON not_in_map.usp=tmp_seed.usp and not_in_map.md5sum=tmp_seed.md5sum) p',
                                   'join taxlite x on p.usp=x.name)'))
  }
  #dbWriteTable(conn,'mdfun_tmp',unique(rseed[,c('md5sum','ufun')]),overwrite=TRUE,temporary=TRUE)
  newUFunmap<-dbGetQuery(conn,paste('select count(*) as cnt from (SELECT distinct md5sum,seedliteid FROM  ( ',
                                    '(SELECT md5sum,x.id as seedliteid FROM tmp_seed join seedlite x on tmp_seed.ufun=x.accession) EXCEPT   ',
                                    '(SELECT d.md5sum as md5sum,seedLiteID FROM seedLite_md5 m join md5 d on m.md5id=d.id) ) not_in_map) mdfun_cnt '))
  if(newUFunmap$cnt>0){
    #dbWriteTable(conn,'funmap_tmp',newUFunmap,overwrite=TRUE,temporary=TRUE)
    rs<-dbSendStatement(conn,paste('insert into seedLite_md5 (md5id,seedLiteID)  (',
                                   'SELECT distinct md5id,seedliteid FROM ',
                                   '( (SELECT md5sum,x.id as seedliteid FROM tmp_seed join seedlite x on tmp_seed.ufun=x.accession join md5 d on  tmp_seed.md5sum=d.md5sum) EXCEPT ',
                                   ' (SELECT md5id,seedLiteID FROM seedLite_md5 m) ) not_in_map)'))
  }
  #dbWriteTable(conn,'abund_tmp',cbind(unique(rseed[,c('md5sum','ab')]),data.frame(smpl=smpl$id)),overwrite=TRUE,temporary=TRUE)
  ab.cnt<-dbGetQuery(conn,'select count(*) as cnt from (select distinct md5sum, mgrastid from tmp_seed except select d.md5sum,s.mgrastid from abundance a join md5 d on a.md5id=d.id join smpl s on a.sampleID=s.id ) ab_cnt')
  if(ab.cnt$cnt>0){
    rs<-dbSendStatement(conn,paste('insert into abundance (md5id,sampleID,abundance) ',
                                   '(select distinct d.id as md5id,s.id as sampleID, ab abundance from tmp_seed t join smpl s on t.mgrastid=s.mgrastid join md5 d on  tmp_seed.md5sum=d.md5sum)'))
  }
}

#' Load  KEGG annotation for particular metagenome to the database
#'
#' @param conn connection to the DB
#' @param mgrastid metagenome id
#' @param path where intermediate files will be stored
#'
#' @export
#'
loadKEGG2DB<-function(conn,mgrastid,path){
  kegg.fname<-normalizePath(gsub('//','/',paste0(path,'/kegg',mgrastid,'.tsv')))
  cat(format(Sys.time(), "%b %d %X"),'loadKEGG',mgrastid,kegg.fname,'\n')
  if(!file.exists(kegg.fname)){
  rkegg<-readKEGG(mgrastid,path)
  cat(format(Sys.time(), "%b %d %X"),'readKEGG',mgrastid,'\n')
  writeDTannot(rannot = rkegg,mgrastID = mgrastid,fname = kegg.fname,annot.col = 'ko')
  }
  beforeTMP<-dbGetQuery(conn,'select count(*) as num from tmp_annot;')
  if(beforeTMP$num>0){
    rs<-dbSendStatement(conn,'delete from tmp_annot;')
  }
  uploadAnnotDT(kegg.fname)
  smpl<-dbGetQuery(conn,'select * from smpl where mgrastid=?',mgrastid)
  first<-dbGetQuery(conn,paste('SELECT count(distinct md5sum) as cnt FROM (',
                               '(SELECT md5sum FROM tmp_annot) EXCEPT (SELECT md5sum FROM asarLite.md5)',
                               ') not_in_md5'))

  newMD5<-first$cnt
  if(newMD5>0){
    newMD5s<-dbGetQuery(conn,paste('SELECT distinct tmp_annot.md5sum FROM (',
                                   '(SELECT md5sum FROM tmp_annot) EXCEPT (SELECT md5sum  FROM asarLite.md5)',
                                   ') not_in_md5 JOIN tmp_annot ON not_in_md5.md5sum=tmp_annot.md5sum'))
    rs<-dbSendStatement(conn,paste('insert into md5 (md5sum) (select distinct md5sum from (',
                                   '(SELECT distinct md5sum FROM tmp_annot) EXCEPT (SELECT md5sum  FROM asarLite.md5)',
                                   ') not_in_md5)'))
  }
  #dbWriteTable(conn,'mdtax_tmp',unique(rseed[,c('md5sum','usp')]),overwrite=TRUE,temporary=TRUE)
  newUSP<-dbGetQuery(conn,paste('SELECT distinct tmp_annot.accession FROM (',
                                '(SELECT accession as access FROM tmp_annot) EXCEPT (SELECT access  FROM asarLite.kort)',
                                ') not_in_tax JOIN tmp_annot ON not_in_tax.access=tmp_annot.accession'))
  if(dim(newUSP)[1]>0){
    write.table(newUSP,'newKEGG.ERROR',append = TRUE)
    cat('WARNING THERE ARE ',dim(newUSP),' UNMAPPED KEGG IDs\n')
  }
  newUSPmap<-dbGetQuery(conn,paste('select count(*) as cnt from (SELECT distinct tmp_annot.md5sum,tmp_annot.accession FROM (',
                                   '(SELECT md5sum,accession FROM tmp_annot) EXCEPT (SELECT d.md5sum as md5sum,k.access as accession FROM md5_kort m ',
                                   'join md5 d on m.md5id=d.id',
                                   'join kort k on m.kortID=k.id)',
                                   ') not_in_map JOIN tmp_annot ON not_in_map.accession=tmp_annot.accession and not_in_map.md5sum=tmp_annot.md5sum) mdtax_cnt '))
  if(newUSPmap$cnt>0){
    rs<-dbSendStatement(conn,paste('insert into md5_kort (md5id,kortID) ',
                                   '(SELECT distinct  p.md5id,x.id as kortid from',
                                   '(SELECT distinct p.md5id,tmp_annot.accession FROM (',
                                   '(SELECT d.id as md5id,md5sum,accession FROM tmp_annot join md5 d on  tmp_annot.md5sum=d.md5sum) EXCEPT (SELECT md5id,d.md5sum as md5sum,k.access as accession FROM md5_kort m ',
                                   'join md5 d on m.md5id=d.id',
                                   'join kort k on m.kortID=k.id)',
                                   ') not_in_map JOIN tmp_annot ON not_in_map.accession=tmp_annot.accession and not_in_map.md5sum=tmp_annot.md5sum) p',
                                   'join kort x on p.accession=x.access)'))
  }
  ab.cnt<-dbGetQuery(conn,'select count(*) as cnt from (select distinct md5sum as md5md5sum, mgrastid from tmp_annot except select md5md5sum,s.mgrastid from abundance a join smpl s on a.sampleID=s.id ) ab_cnt')
  if(ab.cnt$cnt>0){
    rs<-dbSendStatement(conn,paste('insert into abundance (md5id,sampleID,abundance) (',
                                   'select distinct md5sum as md5md5sum,s.id as sampleID, ab as abundance  from ',
                                   '((select distinct md5sum as md5md5sum, mgrastid from tmp_annot) except ',
                                   '(select md5md5sum,s.mgrastid from abundance a join smpl s on a.sampleID=s.id) ) n ',
                                   'join tmp_annot t on t.md5sum=n.md5md5sum and t.mgrastid=n.mgrastid join smpl s on t.mgrastid=s.mgrastid)'))
  }

}

insertStudy2DB<-function(conn,prjID,prjName,prjDesc){
  studyDB<-dbReadTable(conn,'study')
  idx<-which(studyDB$mgrastid==prjID)
  if(length(idx)==0){
      rs<-dbSendStatement(conn,
                          paste('INSERT INTO',
                                'study (name,mgrastid,description)',
                                'VALUES (?,?,?);'),
                          prjName,prjID,prjDesc)
    dbHasCompleted(rs)
    dbGetRowsAffected(rs)
    dbClearResult(rs)
    studyDB<-dbReadTable(conn,'study')
    idx<-which(studyDB$mgrastid==prjID)
  }
  idStudy<-studyDB$id[idx]
  return(idStudy)
}

loadMetagenome<-function(conn,path,mgrastid){
  pDesc<-parseMetadata(paste0(path,mgrastid,'.meta'))
  tbl<-dbGetQuery(conn,'SELECT * FROM study WHERE mgrastid = ?;',
                  pDesc$id)
  if(dim(tbl)[1]==0){
    idStudy<-insertStudy2DB(conn,pDesc$id,pDesc$pname,pDesc$pname)
  }else{
    idStudy<-tbl$id
  }
  insertSample2DB(conn,path,idStudy,mgrastid,pDesc$sname)
}
loadStudy2DB<-function(conn,path,prjID,prjName,prjDesc,updateMV=FALSE){
  mdt<-read.csv(paste0(path,prjID,'.meta.csv'))
  idStudy<-insertStudy2DB(conn,prjID,prjName,prjDesc)
  tbl<-dbGetQuery(conn,'SELECT * FROM mv_smpl_stat WHERE project = ? and cnt>0;',
                  prjID)
  if(dim(tbl)[1]<dim(mdt)[1]){
    if(dim(tbl)[1]==0){
      mdtL<-mdt
    }else{
      idx<-which(!(mdt$MG.RAST.ID%in%tbl$smpl))
      mdtL<-mdt[idx,]
    }
    if(dim(mdtL)[1]>0){
      for(i in 1:dim(mdtL)[1]){
        mgrastid<-as.character(mdtL$MG.RAST.ID[i])
        name<-as.character(mdtL$Metagenome.Name[i])
        insertSample2DB(conn,path,idStudy,mgrastid,name)
        loadSeed2DB(conn,mgrastid,path)
        loadKEGG2DB(conn,mgrastid,path)
      }
      if(updateMV){
        updateMVs(conn)
      }
    }
  }
}


updateMVs <- function(conn){
  # mv_smpl_stat
  rs<-dbSendStatement(conn,'DELETE FROM mv_smpl_stat;')
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Clean mv_smpl_stat\n')
  rs<-dbSendStatement(conn,
                      paste('insert into mv_smpl_stat ',
                            '(stname,project,id,name,smpl,cnt,abundance) ',
                            '(select stname,project,id,name,smpl,cnt,abundance ',
                            'from smpl_stat where cnt>0);'))
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Repopulate mv_smpl_stat\n')

  # mv_dkres
  rs<-dbSendStatement(conn,'DELETE FROM mv_dkres;')
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Clean mv_dkres\n')
  rs<-dbSendStatement(conn,
                      paste('insert into mv_dkres (id,md5,annotation,ko) ',
                            '(select id,md5,annotation,ko from dkres);'))
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Repopulate mv_dkres\n')

  # mv_funtaxab
  rs<-dbSendStatement(conn,'DELETE FROM mv_funtaxab;')
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Clean mv_funtaxab\n')
  rs<-dbSendStatement(conn,
                      paste('insert into mv_funtaxab ',
                            '(md5sum,mgrastid,abundance,usp,ufun) ',
                            '(select md5sum,mgrastid,abundance,usp,ufun from funtaxab);'))
  dbHasCompleted(rs)
  dbGetRowsAffected(rs)
  dbClearResult(rs)
  cat(format(Sys.time(), "%b %d %X"),' Repopulate mv_funtaxab\n')
}
insertSample2DB<-function(conn,path,idStudy,mgrastid,name){
  tbl<-dbGetQuery(conn,'SELECT * FROM smpl WHERE mgrastid = ?',
                  mgrastid)
  if(dim(tbl)[1]==0){
    rs<-dbSendStatement(conn,
                        paste('INSERT INTO',
                              'smpl (name,mgrastid,studyid,ident,length,orgsource,funsource)',
                              'values(?,?,?,?,?,?,?)'),name,mgrastid,idStudy,60,15,"SEED","Subsystems")
    dbHasCompleted(rs)
    dbGetRowsAffected(rs)
    dbClearResult(rs)
  }
}
