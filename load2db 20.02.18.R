#!/usr/bin/Rscript

library(biomformat)
library(RJSONIO)
library(data.table)
library(RCurl)

############ Predefined parameters
 dbname = "asar"
 usr='asar'
 pwd='asar'
 webkey='y8MLy5NRMBfCRGzPcxYP9b86c'
 dtPath = '/Volumes/AS_WD_HFS/OIST/DBData/' #'/var/workspaceR/scalpelData/data/'
 tmpPath = '/Volumes/AS_WD_HFS/Projects/R/MFCMAP/bash/' #'/var/workspaceR/scalpelData/data/'
 rmdPath = '~/Documents/Projects/R/asarDB/load2DB.Rmd'
 rmd.fname = 'load2DB.Rmd'
# archPath = '/var/workspaceR/scalpelData/archive/Burdenko/'
# path<-'/Volumes/AS_WD_HFS/Scalpel/DBData/'
#########################################
source('./path.R')

fl<-dir(path=tmpPath,pattern='project.mgp.*.tbz2$',recursive=TRUE)
wd<-getwd()
for(f in fl){
 d<-sub('.tbz2$','',f)
 if(dir.exists(d)){
   prjTMP<-sub('^project.','',d)
   rm(prjID,prjName,prjDesc)
    server.resource <-
      "http://api.metagenomics.anl.gov/1/project/"
    server.resource <- paste0(server.resource, prjTMP)
    message(paste("Loading the annotations form MG-RAST of",
                  prjTMP),
            domain = NA)
    message("The time spent in this step is proportional to the total amount of remote data...")
    param <-
      list(
        verbosity = 'full',
        auth = webkey
      )
    anno <- tryCatch(
      getForm(
        server.resource,
        .params = param,
        .opts = list(
          noprogress = TRUE)
      ),
      error = function(e) {
        msg <- conditionMessage(e)
        structure(msg, class = "try-error")
      }
    )
    if (inherits(anno, "try-error")) {
      warning(anno)
      return(FALSE)
    }
    invalid.source <- which(grepl("Invalid\\s+ontology\\s+source",
                                  anno))
    if (length(invalid.source))
      stop("invalid ontology source")
    if (length(which(grepl("insufficient\\s+permissions",
                           anno))))
      stop("invalid webkey")
    anno <- fromJSON(
      textConnection(anno),
      header = FALSE,
      sep = "\t",
      stringsAsFactor = F
    )
  prjID<-anno$id
  desc<-as.data.frame(t(as.data.frame(anno$metadata)),stringsAsFactors=FALSE)
  prjName<-desc$project_name
  prjDesc<-desc$project_description
  cat(prjID,prjName,prjDesc,'\n')
  file.copy(rmdPath,paste0(d,'/'))
  setwd(paste0('./',d))
  try(rmarkdown::render(rmd.fname,'pdf_document',clean=FALSE),FALSE,outFile=sub('Rmd$','try.out',rmd.fname))
  setwd(wd)
 }
}

