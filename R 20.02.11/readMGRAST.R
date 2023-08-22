#'Separate functional analysis data into a new table in a file removing duplicates.s
#'Checks md5sum and removes duplicates and copies functional analysis data into the file "d.uspfun".
#'@usage expandNamesDT(x)
#'@param x file that is output by function "mergeAnnots" and from which functional analysis data is going to be extracted.
#'@details Copies coloumns, namely `query sequence id`,`hit m5nr id (md5sum)`,fun and sp with all rows from file "d.merge" to the file "d".
#'@details Changes names of coloumns `query sequence id`,`hit m5nr id (md5sum)`,fun and sp to 'id','md5sum','fun'and 'sp' respectively.
#'@details Checks coloumn 'md5sum' and removes duplicates by transfering rows to "d.ufun" and removes quare brackets.
#'@details Checks coloumn 'md5sum' and 'ufun' and removes duplicates by transfering rows to "d.uspfun" and removes quare brackets.
#'@return table in the file "d.uspfun" which consists of ids, md5sum, function and species name.
#'@seealso @seealso \code{\link{mergeAnnots}}
#'@export
expandNamesDT<-function(d.merge){
  d.merge[,ab:=.N,by=.(`hit m5nr id (md5sum)`)]
  d<-unique(d.merge[,.(`hit m5nr id (md5sum)`,fun,sp,ab)])
  names(d)<-c('md5sum','fun','sp','ab')
  d.uspfun<-d[,expandFunSP.list(fun,sp),by=.(md5sum,ab)]
  return(d.uspfun)
}

expandFunSP.df<-function(fun,sp){
  if(length(fun)>1){
    return(ldply(1:length(fun),.fun=function(.i)expandFunSP.df(fun[.i],sp[.i])))
  }
  ufun<-gsub('(\\]|\\[)','',
             unlist( strsplit( fun, "\\]; *\\[" ) ))
  usp<-gsub('(\\]|\\[)','',unlist( strsplit( sp, "\\]; *\\[" ) ))
  res<-unique(expand.grid(ufun=ufun,usp=usp,fun=fun,sp=sp,stringsAsFactors = FALSE))
  return(res)
}

expandFunSP.list<-function(fun,sp){
  return(as.list(expandFunSP.df(fun,sp)))
}

#'Calculates number of sequences.
#'
#'First calculates sum of sequences ('ab') for every group in bacterial species (usp) and function (ufun) and renames coloumn as 'sum'. Then names 'md5sum as a 'md5' and separates elements of it by comma.
#'Columns 'species' (usp), 'functions' (ufun), 'sum' (sum) and md5 of sequences are returned as a data.table and duplicated rows by all columns are removed.
#'useage getAbundanceMD5FromDT(d.ab)
#'@param d.ab file that should be input for the calculation.
#'@return  file that was input, but adding sum of sequences and sum of md5 in the table.
#'@export
getAbundanceMD5FromDT<-function(d.ab){
  d.ab<-d.ab[,.(sum=sum(ab),md5=paste(md5sum,collapse = ',')),by=.(usp,ufun)]
  d.ab<-unique(d.ab[,.(usp,ufun,sum,md5)])
  return(d.ab)
}

load.data.from.file<-function(path,pattern = '*.3.seed$'){
  path<-sub('/*$','/',path)
  list<-dir(path = path,pattern = pattern)
  cat(paste(list,collapse = '\n'))
  res<-lapply(list,loadMGRAST)
  idx<-which(!sapply(res,is.null))
  res<-res[idx]
  names(res)<-list[idx]
  return(res)
}

#' Checks the presence of GNU head command in the system
#'
#' @return actual command to use as head
#' @export
#'
#' @example
#' makeHeadCMD()
makeHeadCMD<-function(){
  if(Sys.info()["sysname"]=="Darwin"){
    headCMD<-'ghead'
  }else{
    headCMD<-'head'
  }
  system2(headCMD,'--version')->headOut
  if(headOut!=0){
    stop(paste0("Package required command '",headCMD,"' to be installed on the system.\n",
                "MacOS users have to install coreutils by running:\n\t brew install coreutils\n\n"))
  }
  return(headCMD)
}
#' Read single file from MG-RAST folder
#'
#' @param path full path to the file
#' @param fname name of the file
#'
#' @return data.table with content of the file
#' @export
loadMGRAST<-function(path,fname){
  headCMD<-makeHeadCMD()
    l <-readLines(textConnection(system(paste0('tail -n -1 ',path,fname),intern=TRUE)))
  if(!grepl('Download complete',tail(l))){
    warning('Loading of the file "',fname,'"was not complete\n')
    return(NULL)
  }
  res<-fread(cmd=paste0(headCMD,' -n -1 ', paste0(path,fname)),sep='\t',header = TRUE)
  return(res)
}

#'Loading files output by MG-RAST with functional and taxonomical analysis by SEED.
#'
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.kodata.from.file (path)
#' @param path location of a file that should be input.
#' @return list called "ko"
#' @export
#' @example
#' sannot <- load.sdata.from.file()
load.sdata.from.file <- function(path = '.') {
  #  res<-load.data.from.file(path,pattern = '*.3.seed$')
  #  idx<-which(!sapply(res,is.null))
  #  sannot<-res[idx]
  sannot<-load.data.from.file(path,pattern = '*.3.seed$')
  return(sannot)
}

#'Loading files output by MG-RAST with functional analysis by SEED.
#'
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.fdata.from.file(path)
#' @param path location of a file that should be input.
#' @return list called "fannot"
#' @export
#' @example
#' fannot <-load.fdata.from.file(path = ".")
load.fdata.from.file <- function(path = ".") {
  fannot<-load.data.from.file(path,pattern = '*.3.fsub$')
  return(fannot)
}

#'Loading file KEGG Orthology output by MG-RAST.
#'
#' After entering location of a file (path) and name of the file (pattern), this function reads in information from the file as a list.
#' @usage load.kodata.from.file (path)
#' @param path location of a file that should be input.
#' @return list called "ko"
#' @export
#' @example
#' ko <-load.kodata.from.file()
load.kodata.from.file <- function(path = '.') {
  ko<-load.data.from.file(path,pattern = '^m.*.ko$')
  return(ko)
}


make.d.res <- function(.list) {
  d.res<-ldply(.data = .list$res,
               .fun = function(.x){
                 ab<-.x$ab;
                 ab$mgid=.x$name;
                 return(ab)
               })
}

#command "d.kres <- make.d.kres(kres.res)" should be run
make.d.kres <- function(.list){
  d.kres<-unique(ldply(.data = .list$kres,
                       .fun = function(.x){
                         ab<-.x$ab;
                         return(ab)
                       }))
}

our.aggregate <- function(d.res) {
  data.table::dcast(setDT(d.res),usp+ufun+md5 ~ mgid,value.var = 'sum',fill = 0)->d.bm

  d.sp<-aggregate(.~usp,as.data.frame(d.bm)[,-c(2,3)],FUN = sum)
  d.fun<-aggregate(.~ufun,as.data.frame(d.bm)[,-c(1,3)],FUN = sum)
  return(d.bm)
}

#command "kres.res <- our.merge()" should be run
our.merge <- function(path='.') {
  res<-list()
  kres<-list()
  fannot <- load.fdata.from.file(path = path)
  sannot <- load.sdata.from.file(path = path)
  ko <- load.kodata.from.file(path = path)
  nms<-gsub('.fsub$','',names(fannot))
  for(i in 1:length(fannot)){
    f<-fannot[[i]]
    f$`query sequence id`<-gsub('\\|Subsystems$','',f$`query sequence id`)
    s<-sannot[[i]]
    s$`query sequence id`<-gsub('\\|SEED$','',s$`query sequence id`)
    d.k<-unique(ko[[i]][,list(`hit m5nr id (md5sum)`,
                              `semicolon separated list of annotations`)])
    #creates 3rd column with accession number itself only
    d.k1<-d.k[,list(`semicolon separated list of annotations`,
                    ko=unlist(gsub('accession=\\[K([0-9]+)\\].*','K\\1',
                                   unlist(str_split(`semicolon separated list of annotations`,
                                                    ';'))))),
              by=.(`hit m5nr id (md5sum)`)]
    names(d.k1)<-c('md5','annotation','ko')
    #kres is ready
    kres[[nms[i]]]<-list(ab=d.k1,name=nms[i])

    keycols<-names(f)[1:12]
    setkeyv(f,keycols)
    setkeyv(s,keycols)
    d.merge<-merge(f,s,all=TRUE,suffixes = c('.sub','.sp'))
    names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
    d.m1<-d.merge[,list(sub,sp,`query sequence id`,fun=unlist(gsub('accession=\\[SS([0-9]+)\\].*','SS\\1',
                                                                   unlist(str_split(sub,';'))))),
                  by=.(`hit m5nr id (md5sum)`)]
    d.uspfun<-expandNamesDT(d.m1)
    d.ab<-getAbundanceMD5FromDT(d.uspfun)
    #res is ready
    res[[nms[i]]]<-list(ab=d.uspfun,name=nms[i])

    #d.m.t<-table(d.merge$fun,d.merge$sp)
    cat(paste(i,nms[i],'\n'))
  }
  return(list(kres=kres, res=res))
}

#' Read SEED files and prepare aggregated data.dable
#'
#' @param mgrastid metagenome ID to deal with
#' @param path path to the download folder
#'
#' @return data.table with all reads aggregated at the MD5 level
#' @export
readSEED<-function(mgrastid,path){
  path<-sub('/*$','/',path)
  s<-loadMGRAST(path,paste0(mgrastid,'.seed'))
  sa<-unique(s[,.(ab=length(`query sequence id`),`semicolon separated list of annotations`),by=.(`hit m5nr id (md5sum)`)])
  f<-loadMGRAST(path,paste0(mgrastid,'.fsub'))
  fa<-unique(f[,.(ab=length(`query sequence id`),`semicolon separated list of annotations`),by=.(`hit m5nr id (md5sum)`)])
 # f$`query sequence id`<-gsub('\\|Subsystems$','',f$`query sequence id`)
#  s$`query sequence id`<-gsub('\\|SEED$','',s$`query sequence id`)
  keycols<-names(fa)[1:2]
  setkeyv(fa,keycols)
  setkeyv(sa,keycols)
  d.merge<-merge(f,s,all=TRUE,suffixes = c('.sub','.sp'))
  names(d.merge)<-gsub('semicolon separated list of annotations.','',names(d.merge))
  d.m1<-d.merge[,list(sub,sp,ab,fun=unlist(gsub('accession=\\[SS([0-9]+)\\].*','SS\\1',
                                                                 unlist(str_split(sub,';'))))),
                by=.(`hit m5nr id (md5sum)`)]
  d.uspfun<-d.m1[,expandFunSP.list(fun,sp),by=.(`hit m5nr id (md5sum)`,ab)]
  names(d.uspfun)<-c("md5sum","ab","ufun","usp","fun","sp")
  return(d.uspfun)
}

readKEGG<-function(mgrastid,path){
  path<-sub('/*$','/',path)
  ko <- loadMGRAST(path,paste0(mgrastid,'.ko'))
  ko[,ab:=.N,by=.(`hit m5nr id (md5sum)`)]
  d.k<-unique(ko[,list(`hit m5nr id (md5sum)`,ab,
                            `semicolon separated list of annotations`)])
  #creates 3rd column with accession number itself only
  d.k1<-d.k[,list(`semicolon separated list of annotations`,ab,
                  ko=unlist(gsub('accession=\\[K([0-9]+)\\].*','K\\1',
                                 unlist(str_split(`semicolon separated list of annotations`,
                                                  ';'))))),
            by=.(`hit m5nr id (md5sum)`)]
  d.k2<-d.k1[grep('^K[0-9]+$',ko)]
  names(d.k2)<-c('md5sum','annotation','ab','ko')
  return(d.k2)
}

#' Title
#'
#' @param prjTMP permanent (mgp) or temporary MG-RAST project ID
#' @param webkey valid web-key from MG-RAST
#'
#' @return
#' @export
#'
#'
loadProjectCSV<-function(prjTMP,webkey){
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
  proj.ID<-anno$id
  metagenomes <-
    ldply(anno$metagenomes, function(.x) {
      data.frame(
        MG.RAST.ID = .x$metagenome_id,
        Metagenome.Name = .x$name,
        bp.Count = .x$basepairs,
        Sequence.Count = .x$sequences,
        Biome = .x$biome,
        Feature = .x$feature,
        Material = .x$material,
        Location = .x$location,
        Country = .x$country,
        Sequence.Type = .x$sequence_type,
        Sequence.Method = .x$sequencing_method
      )
    })
  write.csv(metagenomes,file = paste0(proj.ID,'.meta.csv'))
  return(metagenomes)
}

parseMetadata<-function(meta.file){
  p<-fromJSON(meta.file)
  return(list(id=p$project$id,pname=p$project$name,sname=p$sample$name))
}