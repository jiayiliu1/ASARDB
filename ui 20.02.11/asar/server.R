library(shiny)
library(ggplot2)
#library(mmnet)
library(gplots)
library(data.table)
library(plyr)
library(dplyr)
library(pathview)
library(stringr)
library(biomformat)
library(KEGGREST)
library(DT)
library(shinythemes)
library(matrixStats)
library(png)  # For writePNG function
#library(devtools)
library(d3heatmap)
library(edgeR)
library(RColorBrewer)
library(rhandsontable)

library(limma)
#source('mg.key.R')
load("pathview.Rdata")
ls()  # List objects in the current environment

#load('~/Documents/Projects/R/MFCMAP/R/pathview.Rdata')

library(asarDB)

options(shiny.maxRequestSize=10*1024^3)
#mdt$Group<-mdt$Metagenome.Name

ko.path.name<-ko.path.name[-grep('ko01100',ko.path.name$ko),]
kegg <- kegg[-grep('ko01100', kegg$ko),]


#' Title
#'
#' @param result2 -- data matrix with only interesting subset of samples
#' @param t1  -- taxonomy selection level
#' @param tn  -- taxonomy selection value
#' @param f1  -- function selection level
#' @param fn  -- function selection level
#' @param t2  -- taxonomy aggregation level
#' @param f2  -- function aggregation level
#' @param fun  -- aggregation function (default sum)
#'
#' @return
#' @export
#'
#' @examples
Intfuntax <- function(result2, t1, tn, f1, fn, t2=NULL, f2=NULL,s=NULL,fun=sum,aggSamples=c()){
  cat('t1=',t1, 'tn=',paste0('|',tn,'|'), 'f1=',f1, 'fn=',fn, 't2=',t2, 'f2=',f2,'\n')
  cat(dim(result2),'\n',names(result2),'\n')
  #selagg <- function(data, taxS, taxVal, funS, funVal,samples, taxAgg=NULL, funAgg=NULL,aggSamples=c(),aggFunction=sum){
  result2<-selagg(result2,taxS=t1,taxVal=tn,funS=f1,funVal=fn,taxAgg=t2,funAgg=f2,samples=s,aggSamples=aggSamples)
  # if(t1!="toplevel"){
  #   #quo_tn <- enquo(t1)
  #   # idx<-unique(unlist(sapply(tn,grep,x=result2[,get(t1)])))
  #   # result2 <- filter(result2,grepl(tn,!!quo_tn))
  #   #result2<-dplyr::filter_(result2,paste0(t1,"=='",tn,"'"))
  #   result2<-dplyr::slice_(result2,paste0("grep('",tn,"',",t1,')'))
  # }
  # if(f1!="toplevel"){
  #   #result2 <- result2[grep(fn, result2[,get(f1)])]
  #   result2<-dplyr::slice_(result2,paste0("grep('",fn,"',",f1,')'))
  # }
  # if(!is.null(t2) & !is.null(f2)){
  #   result2<- plyr::ddply(result2, c(t2,f2), numcolwise(fun))
  # }else{
  #   if(!is.null(t2)){
  #     result2<- plyr::ddply(result2, t2, numcolwise(fun))
  #   }
  #   if(!is.null(f2)){
  #     result2<- plyr::ddply(result2, f2, numcolwise(fun))
  #   }
  # }
  # cat(dim(result2),'\n',names(result2),'\n')
  #
  # result2 <- data.table(result2[which(!is.na(result2[,1])& !is.na(result2[,2])& result2[,2]!= ""),])
  # cat(dim(result2),'\n',names(result2),'\n')
  return(result2)
}
# Normalization function using TMM
normalize_TMM <- function(obj) {
  dge <- DGEList(counts = obj)
  dge <- calcNormFactors(dge)
  obj_norm <- cpm(dge, log = FALSE)
  return(obj_norm)
}

# Normalization function using DESeq2's median of ratios
normalize_DESeq2 <- function(obj) {
  # Compute the row medians
  row_medians <- apply(obj, 1, median, na.rm = TRUE)

  # Check for zero row medians
  zero_medians <- row_medians == 0
  row_medians[zero_medians] <- 1 # This is just a placeholder to prevent divide by zero

  # Calculate median ratio
  med_ratio <- median(row_medians, na.rm = TRUE) / row_medians

  # Check for NA in med_ratio
  if(any(is.na(med_ratio))) {
    stop("NA values encountered in median ratios")
  }

  # Normalize the object
  obj_norm <- t(t(obj) * med_ratio)

  return(obj_norm)
}


# Normalization function using RPKM/FPKM method
normalize_RPKM <- function(obj) {
  # Compute total reads
  total_reads <- rowSums(obj, na.rm = TRUE)

  # Compute effective lengths
  effective_lengths <- colSums(obj > 0, na.rm = TRUE)

  # Check for zero effective lengths (and potentially replace or warn)
  zero_lengths <- effective_lengths == 0
  if(any(zero_lengths)) {
    warning("Some samples have zero effective lengths. Setting to 1 to avoid division by zero.")
    effective_lengths[zero_lengths] <- 1 # This is just a placeholder to prevent divide by zero
  }

  # Normalize the object
  obj_norm <- 1e9 * t(t(obj) / (total_reads / 1e6 * effective_lengths))

  # Check for NA/NaN/Inf values and handle or warn
  if(any(is.na(obj_norm) | is.nan(obj_norm) | is.infinite(obj_norm))) {
    warning("NA/NaN/Inf values detected after normalization.")
  }

  return(obj_norm)
}



# # Normalization function using DESeq2's median of ratios
# normalize_DESeq2 <- function(obj) {
#   row_medians <- apply(obj, 1, median, na.rm = TRUE)
#   med_ratio <- median(row_medians, na.rm = TRUE) / row_medians
#   obj_norm <- t(t(obj) * med_ratio)
#   return(list(data=obj_norm, factors=med_ratio))
# }
#
# # Normalization function using RPKM/FPKM method
# normalize_RPKM <- function(obj) {
#   total_reads <- rowSums(obj, na.rm = TRUE)
#   effective_lengths <- colSums(obj > 0, na.rm = TRUE)
#   rpkm_factors <- total_reads / (1e6 * effective_lengths)
#   obj_norm <- 1e9 * t(t(obj) / rpkm_factors)
#   return(list(data=obj_norm, factors=rpkm_factors))
# }
#
# # Normalization function using TMM
# normalize_TMM <- function(obj) {
#   dge <- DGEList(counts = obj)
#   dge <- calcNormFactors(dge)
#   obj_norm <- cpm(dge, log = FALSE)
#   return(list(data=obj_norm, factors=dge$samples$norm.factors))
# }



# # Normalization function using DESeq2's median of ratios
# normalize_DESeq2 <- function(obj) {
#   row_medians <- apply(obj, 1, median, na.rm = TRUE)
#   med_ratio <- median(row_medians, na.rm = TRUE) / row_medians
#   obj_norm <- t(t(obj) * med_ratio)
#   return(obj_norm)
# }
#
# # Normalization function using RPKM/FPKM method
# normalize_RPKM <- function(obj) {
#   total_reads <- rowSums(obj, na.rm = TRUE)
#   effective_lengths <- colSums(obj > 0, na.rm = TRUE)
#   obj_norm <- 1e9 * t(t(obj) / (total_reads / 1e6 * effective_lengths))
#   return(obj_norm)
# }

# Plot density of normalized values
plot_density <- function(normalized_data, method_name) {
  require(ggplot2)
  data_long <- as.data.frame(as.matrix(normalized_data))
  # Convert data to long format for ggplot
  data_long <- reshape2::melt(normalized_data)

  if(ncol(data_long)==1){
    p <- ggplot(data_long, aes(x=value)) +
      geom_density(fill="blue", alpha=0.5) +
      labs(title=paste("Density plot for", method_name, "normalized data"),
           x="Normalized values",
           y="Density")
  }else{
    p <- ggplot(data_long, aes(x=value, fill=variable)) +
      geom_density(alpha=0.5) +
      labs(title=paste("Density plot for", method_name, "normalized data"),
           x="Normalized values",
           y="Density")
  }

  return(p)
}

# Calculate summary statistics for normalized data
summary_stats <- function(normalized_data) {
  if (is.null(ncol(normalized_data)) == TRUE) {
    mean_val = mean(normalized_data, na.rm=TRUE)
    median_val = median(normalized_data, na.rm=TRUE)
    sd_val = sd(normalized_data, na.rm=TRUE)
    cv_val = sd_val / mean_val

    return(list(mean=round(mean_val,3), median=round(median_val,3), sd=round(sd_val,3), cv=round(cv_val,3)))
  } else if (ncol(normalized_data) == 2) {
    mean_val1 = mean(normalized_data[,1], na.rm=TRUE)
    mean_val2 = mean(normalized_data[,2], na.rm=TRUE)
    median_val1 = median(normalized_data[,1], na.rm=TRUE)
    median_val2 = median(normalized_data[,2], na.rm=TRUE)
    sd_val1 = sd(normalized_data[,1], na.rm=TRUE)
    sd_val2 = sd(normalized_data[,2], na.rm=TRUE)
    cv_val1 = sd_val1 / mean_val1
    cv_val2 = sd_val2 / mean_val2

    return(list(mean_col1=round(mean_val1,3), mean_col2=round(mean_val2,3),
                median_col1=round(median_val1,3), median_col2=round(median_val2,3),
                sd_col1=round(sd_val1,3), sd_col2=round(sd_val2,3),
                cv_col1=round(cv_val1,3), cv_col2=round(cv_val2,3)))
  } else {
    stop("Input data has more than two columns")
  }
}


loadRdata <- function(fname) {
  if(file.exists(fname)){
    load(file = fname)
    #!TODO: Add misnumber check and report to the user
    mdt <- mdt[which(!is.na(match(rownames(mdt), colnames(funtaxall)))), ]
    idx <- match(getSampleNames(funtaxall), rownames(mdt))
    cat('loadData-1:',idx,'\n')
    mdt<-mdt[idx,]
    cat('loadData-2:\t',rownames(mdt),'\n\t\t',getSampleNames(funtaxall),'\n')
    # Apply normalization to the data if the method is specified and not "none"
    if (normalize() != "none") {
      norm_method <- normalize
      mdt <- returnAppropriateObj(mdt, norm = TRUE, log = TRUE, method = norm_method)
    }

    list <-
      list(
        "funtaxall" = funtaxall,
        "mdt" = mdt,
        "kegg" = kegg,
        "ko.path.name" = ko.path.name,
        "d.kres" = d.kres
      )
  }else{
    list<-list()
  }
  return(list)
}


getFunction<-function(funStr){
  cat('getFunction:',funStr,'\n')
  switch (funStr,
          'sum' = sum,
          'sd' = sd,
          'min' = min,
          'max' = max,
          sd
  )
}
make2d <- function(funtax,numRows=50,sub.na=0){
  obj <- matrix(nrow = length(unique(funtax[,2])), ncol = length(unique(funtax[,1])))
  colnames(obj) <- unique(funtax[,1])
  rownames(obj) <- unique(funtax[,2])
  for (x in 1:nrow(funtax)){
    obj[as.character(funtax[x,2]), as.character(funtax[x,1])]<- funtax[x,3]
  }
  obj[is.na(obj)]<-sub.na
  rowmean <- rowMeans(obj)
  colmean <- colMeans(obj)
  obj <- obj[which(rowmean$Means != 0), which(colmean$Means != 0)]
  return(obj)
}

plotHeatmap <- function(obj, n, normMethod = "none", log = TRUE, fun = sd, method = 'TMM', ...) {
  mat = returnAppropriateObj(obj, norm = (normMethod != "none"), log, method)
  cat(dim(mat),'\n')
  otusToKeep = which(rowSums(mat) > 0)
  cat(otusToKeep,'\n')
  if(length(otusToKeep)==0){
    return(matrix(-1,ncol = 1,nrow = 1))
  }else if(length(otusToKeep)==1){
    otuStats<-fun(mat)
  }else{
    otuStats = apply(mat[otusToKeep, ], 1, fun)
  }
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(length(otusToKeep),n,dim(mat)[1]))]]
  mat2 = mat[otuIndices, ]
  return(mat2)
}

makePlotTitle<-function(title,tl1, tn, fl1, fn, tl2=NULL, fl2=NULL,s=NULL,s2='MGRAST'){
  if(fl1=='toplevel') fn<-'root'
  if(tl1=='toplevel') tn<-'root'
  main1 <-
    paste(
      title,'\n',
      'Metagenome dimension:\n',
      '\tSelection:\t',paste0(s,collapse = ', '),'\n',
      '\tAggregation:\t',s2,'\n')
  main2 <-
    paste('Functional dimension:\n',
          '\tSelection:\n\t\tlevel:\t"',fl1,'"\n\t\tvalues:\n',
          paste0('\t\t\t"', fn, '"\n'))
  main3 <- ifelse(is.null(fl2),'',
                  paste('\tAggregation:\n\t\tlevel:\t"',fl2,'"\n'))
  main4 <-
    paste('Taxonomy dimension:\n',
          '\tSelection:\n\t\tlevel:\t"',tl1,'"\n\t\tvalues:\n',
          paste('\t\t\t"',tn,'"\n', collapse = ', '))
  main5 <- ifelse(is.null(tl2),'',
                  paste('\tAggregation:\n\t\tlevel:\t"',tl2,'"\n'
                  ))
  return(paste(main1,main2,main3,main4,main5))
}

returnAppropriateObj <- function(obj, norm, log, method) {
  obj1=as.data.frame(obj)
  obj=as.matrix(obj1[,unlist(lapply(obj1,is.numeric))])
  if (class(obj)[1] != 'matrix' | any(dim(obj) == 0)) {
    return(matrix(NA, ncol = 1, nrow = 1))
  }

  if (log) {
    obj <- log2(obj + 1)
  }

  if (norm) {
    switch(method,
           'TMM' = {
             obj <- normalize_TMM(obj)
           },
           'DESeq2' = {
             obj <- normalize_DESeq2(obj)
           },
           'RPKM' = {
             obj <- normalize_RPKM(obj)
           }
    )
  }
  d=data.frame(cbind(obj1[,!unlist(lapply(obj1,is.numeric))],obj))
  colnames(d)=colnames(obj1)

  return(data.table(d))
}
# returnAppropriateObj <- function(obj, norm, log, method) {
#   obj1 = as.data.frame(obj)
#   obj = as.matrix(obj1[, unlist(lapply(obj1, is.numeric))])
#
#   if (class(obj)[1] != 'matrix' | any(dim(obj) == 0)) {
#     return(matrix(NA, ncol = 1, nrow = 1))
#   }
#
#   if (log) {
#     obj <- log2(obj + 1)
#   }
#
#   norm_factors <- NULL  # Placeholder for normalization factors
#
#   if (norm) {
#     switch(method,
#            'TMM' = {
#              norm_result <- normalize_TMM(obj)
#              obj <- norm_result$data
#              norm_factors <- norm_result$factors
#            },
#            'DESeq2' = {
#              norm_result <- normalize_DESeq2(obj)
#              obj <- norm_result$data
#              norm_factors <- norm_result$factors
#            },
#            'RPKM' = {
#              norm_result <- normalize_RPKM(obj)
#              obj <- norm_result$data
#              norm_factors <- norm_result$factors
#            }
#     )
#   }
#
#   d = data.frame(cbind(obj1[, !unlist(lapply(obj1, is.numeric))], obj))
#   colnames(d) = colnames(obj1)
#
#   # Returning both the normalized data and the normalization factors
#   return(list(data = data.table(d), factors = norm_factors))
# }

############### KEGG Utility functions #########
getpathsfromKOs <- function(KOs){
  setkey(kegg, K)
  temp <- kegg[KOs]
  unlist(str_split(gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ",")), ','))
}
getPathwayList <- function(funtax, sp.li, mgm, ko_sd) {
  dk7 <- filter_stats(funtax, sp.li, mgm, ko_sd)
  kos<- unique(dk7[,"ko"])
  listko <- getpathsfromKOs(unique(dk7[,"ko"]))
  ko.path.name$ko <- gsub('ko','',ko.path.name$ko)
  pathandnames <- as.matrix(ko.path.name[which(ko.path.name$ko %in% listko),])
  return(pathandnames)
}
getpathfromKO <- function(KO){
  cat(unlist(KO),'\n')
  temp <- kegg[K == KO]
  cat(dim(temp),'\n')
  temp <- gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ","))
  cat(length(temp),'\n')
  return(temp)
}
get_ko_data <- function(funtax, taxon, metagenomes) {
  idx <- unique(unlist(sapply(taxon, grep, x = funtax$usp)))
  if(length(idx)>0){
    d<-funtax[idx]#getSpecieFromAbundMD5_2(funtax,sp = taxon,aggregate = FALSE)
  }else{
    d<-funtax
  }
  cat('get_ko_data-1:',dim(d),names(d),'\n')
  indC<-c(which(names(d)=='md5'),match(metagenomes,names(d)))
  cat('get_ko_data-2:',length(indC),indC,'\n')
  #  d5.1<-d[,list(m5=unlist(str_split(md5,','))),by=.(usp,ufun,md5)]
  d5.1<-unique(d[,list(m5=unlist(str_split(md5,','))),by=.(ufun,md5)])
  cat('get_ko_data-3:',dim(d5.1),names(d5.1),'\n')
  d5<-merge(d5.1,d[,..indC],by='md5')
  cat('get_ko_data-4:',dim(d5),names(d5),'\n')
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  cat('get_ko_data-5:',dim(dk5),names(dk5),'\n')
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'ufun', 'annotation')]),FUN=sum)
  cat('get_ko_data-6:',dim(adk5),names(adk5),'\n')
  return(adk5)
}

pathImage<-function(funtax, sp.li, mgm, pathwi, kostat,nms) {
  withProgress(message = paste("Drawing KEGG pathway", pathwi, "for", sp.li, ".", "Please wait!"), value = 10, {
    adk5<-get_ko_data(funtax, sp.li, mgm)
    rownames(adk5)<-adk5$ko
    adk5<-adk5[,-1]
    indRow <- match(rownames(adk5),rownames(kostat))
    indCol <- match(colnames(adk5),colnames(kostat))
    if(any(is.na(indCol))|any(is.na(indRow))){
      showModal(modalDialog(
        title = titleForNonMatchKostat, textForNonMatchKostat, easyClose = TRUE, footer = NULL
      ))
    } else {
      kostat <- kostat[indRow, indCol]
      ind0 <- apply(kostat,2,function(.x){which(.x==0)})
      if(any(sapply(1:length(ind0), function(i){any(adk5[ind0[[i]],names(ind0)[i]]!=0)}))){
        showModal(modalDialog(
          title = titleForNonMatchZerosKostat, textForNonMatchZerosKostat, easyClose = TRUE, footer = NULL
        ))
      } else {
        sapply(1:length(ind0), function(i){kostat[ind0[[i]],names(ind0)[i]] <<- 10^(-5)})
        adk5 <- adk5/kostat*100
        obj<-as.matrix(adk5)
        colnames(obj)<-nms
        obj<-avearrays(obj)
        idx<-match(kegg$K[kegg$ko==paste0('ko',pathwi)],rownames(obj))
        idx<-idx[!is.na(idx)]
        if(length(idx)>0&
           diff(range(as.vector(obj[idx,])))
        ){
          save(obj,pathwi,sp.li,file=paste0('dump.',pathwi,'.',sp.li,'.Rdata'))
          pathview(gene.data = obj, pathway.id = pathwi,
                   species = "ko", out.suffix = paste0(sp.li,".ko"), kegg.native = T,
                   limit = list(gene=range(as.vector(obj[idx,])),cpd=1))
          mg.key(fname = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png"),
                 names = colnames(obj),node.size = node.size,cex = 0.25,lwd=0.25,tstep = 0.05,rstep = 0.05)
        }else{
          showModal(modalDialog(
            title = titleForNonMatchZerosKostat, textForNonMatchZerosKostat, easyClose = TRUE, footer = NULL
          ))
        }
      }
    }
  })
}
filter_stats <- function(funtax, taxon, metagenomes, sd_cutoff) {
  adk5 <- get_ko_data(funtax, taxon, metagenomes)
  cat('filter_stats-1:',dim(adk5),names(adk5),'\n')
  indM <- which(names(adk5)%in%c(metagenomes))
  cat('filter_stats-2:',length(indM),indM,'\n')
  adk5 <- as.data.table(adk5)
  cat('filter_stats-3:',taxon, metagenomes, sd_cutoff,'\n')
  dk6 <- data.frame(ID = adk5[,"ko"], Means=rowMeans(adk5[,..indM]), SD=rowSds(as.matrix(adk5[,..indM])))
  cat('filter_stats-4:',dim(dk6),names(dk6),'\n')
  dk6idx<-order(dk6$SD,decreasing = TRUE)[which((dk6$Means!=0))]
  cat('filter_stats-5:',length(dk6idx),'\n')
  cutoff<-as.integer(max(2,sd_cutoff*length(dk6idx)/100))
  cat('filter_stats-6:',cutoff,'\n')
  dk7 <- adk5[dk6idx[1:cutoff],]
  cat('filter_stats-7:',dim(dk7),names(dk7),'\n')
  return(dk7)
}
# getSpecieFromAbundMD5<-function(taxall, tx=tx, sp = SpName, aggregate=FALSE){
#   drops <- c("domain","phylum", "class", "order", "family", "genus", "species" ,"usp" )
#   drops <- drops[drops!= tx]
#   dd.res <- taxall[ , !(names(taxall) %in% drops), with = FALSE]
#   idx<-unique(unlist(sapply(sp,grep,x=dd.res[,get(tx)])))
#   d.res<-dd.res[idx]
#   if(aggregate&dim(d.res)[1]>1) d.res<-aggregate(.~ufun,as.data.frame(d.res)[,-c(1,2)],FUN = sum)
#   return(d.res)
# }
# getSpecieFromAbundMD5_2<-function(d.bm,sp=SpName,aggregate=FALSE){
#   idx<-unique(unlist(sapply(sp,grep,x=d.bm$usp)))
#   d.res<-d.bm[idx]
#   if(aggregate&dim(d.res)[1]>1) {
#     d.res<-aggregate(.~ufun,
#                      as.data.frame(d.res)[,-c(1,3)],
#                      FUN = sum)
#     }
#   return(d.res)
# }
# getpathfromKO <- function(KO){
#   temp <- kegg[K == KO]
#   temp <- gsub('ko','',paste0(unlist(temp[,"ko"]), collapse = ","))
# }
pathwayHeatmap<-function(funtax,sp.lis, mgms, ko_sd) {
  adk5 <- filter_stats(funtax, sp.lis, mgms, ko_sd)
  cat(dim(adk5),'\n')
  lastcol<- ncol(adk5)+1
  for (y in 1:nrow(adk5)){
    adk5[y,"pathwayID"] <- getpathfromKO(adk5[y,"ko"])
  }
  indC<-which(names(adk5)%in%c('ko', mgms))
  adk5 <- as.data.table(adk5)
  adk6<-adk5[,list(pat=unlist(str_split(pathwayID,','))),by=.(ko)]
  a6<-merge(adk6,adk5[,..indC],by='ko')
  a7<-aggregate(.~pat,as.data.frame(a6[,-c('ko')]),FUN=sum)
  pnameIdx<-match(paste0('ko',a7$pat),ko.path.name$ko)
  a7<-a7[!is.na(pnameIdx),]
  pnameIdx<-pnameIdx[!is.na(pnameIdx)]
  rownames(a7)<-ko.path.name$name[pnameIdx]
  a8<- as.matrix(a7[,-1])
  rownames(a8) <- ko.path.name$name[pnameIdx]
  return(a8)
}

kostat <- function(funtaxall, d.kres){
  funtaxall <- funtaxall[,-c("usp","species", "genus", "family", "order", "class", "phylum", "domain","ufun", "FUN2", "FUN3", "FUN4")]
  d5.1<-funtaxall[,list(m5=unlist(str_split(md5,','))),by=.(md5)]
  d5<-unique(merge(d5.1,funtaxall,by='md5'))
  dk5<-unique(merge(d5,d.kres,all=FALSE,by.x='m5',by.y='md5')[,-c('md5','.id')])
  adk5<-aggregate(.~ko,as.data.frame(dk5[,-c('m5', 'annotation')]),FUN=sum)
  rownames(adk5)<- adk5$ko
  adk5 <- adk5[,-1]
}


# Define server logic required to draw a histogram
shinyServer(function(input, output, session) {

  ############### Utility functions #########
  chooseDends <- function(res){
    if(dim(res)[1]==1&dim(res)[2]==1){
      dend <- 'none'
    } else if(dim(res)[1]==1){
      dend <- 'column'
    } else if(dim(res)[2]==1){
      dend <- 'row'
    } else{
      dend <- 'both'
    }
    return(dend)
  }

  ############### Reactive expressions #########
  set_taxlevel1 <-reactive({input$set_taxlevel1})
  set_taxlevel2 <-reactive({input$set_taxlevel2})
  set_funlevel1 <-reactive({input$set_funlevel1})
  set_funlevel2 <-reactive({input$set_funlevel2})
  colorPalette <-reactive({input$colorPalette})
  colPal <- reactive({brewer.pal(9, colorPalette())})
  colName <- reactive({input$colName})
  colNameS <- reactive({colnames(mdt)})
  sampleNames<-reactiveVal(value=setNames(getSampleNames(funtaxall), mdt[,1]),label='sampleNames')
  mgall <-reactive({input$mgall})   # list of samples of interest
  mg1   <-reactive({input$mg1})     # one sample of interest for F/T tab
  mg2   <-reactive({input$mg2})     # one sample of interest for F/M tab
  mg3   <-reactive({input$mg3})     # one sample of interest for T/M tab
  tl1   <-reactive({input$tl1})     # taxonomy selectinon level
  tl2   <-reactive({input$tl2})     # taxonomy aggregation level
  fl1   <-reactive({input$fl1})     # function selection level
  fl2   <-reactive({input$fl2})     # function aggregation level
  trf   <-reactive({input$trfun})     # function for heatmap truncation
  trnum   <-reactive({
    if(grepl('^[0-9.,]+$',input$trnum)){
      as.integer(input$trnum)
    }else{
      showModal(
        modalDialog(
          title = titleForNRowErrorPopup,
          textForNRowErrorPopup,
          easyClose = TRUE,
          footer = NULL
        )
      )
      50
    }
  })     # number of rows for heatmap truncation
  taxnames<-reactiveValues(
    tn=as.character(funtaxall$genus[(nrow(funtaxall)/2)]),
    fn=as.character(funtaxall$FUN4[(nrow(funtaxall)/2)]))

  tn<-reactive({                   #
    if(tl1()=='toplevel'){
      taxnames$tn<-'toplevel'
    }else{
      taxnames$tn<-input$tn
    }
    return(taxnames$tn)
  })

  pathw <- reactive({input$PathwayID})
  kostat <- kostat(funtaxall, d.kres)
  fn    <-reactive({input$fn})     #
  ko_sd <-reactive({input$ko_sd})  #
  numrow1 <- reactiveValues(plot1 =0,plot2 =0,plot3 =0,plot4 =0)
  downHeat1 <- reactiveValues(is =TRUE)
  downHeat2 <- reactiveValues(is =TRUE)
  downHeat3 <- reactiveValues(is =TRUE)
  downHeat4 <- reactiveValues(is =TRUE)
  plotRes <- reactiveValues(plot1=NULL,plot2=NULL)

  # Add the normalization method selection input
  norm_method <- reactive({
    input$normMethod
  })

  ############### Heatmap Plot1 #########
  output$plot1 <- renderD3heatmap({
    cat('output$plot1\n')
    tl1 <- tl1()  # taxonomy selection level
    tl2 <- tl2()  # taxonomy aggregation level
    tn  <- tn()   # taxonomy selection value
    fl1 <- fl1()  # function selection level
    fl2 <- fl2()  # function aggregation level
    fn  <- fn()   # function selection value
    mg1 <- mg1()  # selected sample value
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation
    cat('plot1-0.1:',trf(),trnum,input$trfun,mg1,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat1$is <- FALSE
      return()
    }
    #  cat(match(fl1,dimNames('tax')),match(fl2,dimNames('tax')))
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))
    cat(paste0('|',c(tl1, tl2, fl1, fl2, mg1),'|'),'\n')#,paste0('|',names(funtaxall),'|'),'\n')
    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
    #    funtax <- funtaxall[, ..keepcols]
    #    funtax <- dplyr::select(funtaxall,c(tl1, tl2, fl1, fl2, mg1))
    #    cat(names(funtax),'\n',tl1, tn, fl1, fn, tl2, fl2,'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2,s=mg1))



    # Apply the selected normalization method to each matrix in the list

    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
      # normalization_factors <- funtax
      # print(normalization_factors)
      #print(333)
    }
    funtax<<-funtax
    cat(input$fl2)
    cat(input$tl2)
    cat(input$trnum)
    cat(12345)
    trfun<<-getFunction(trf())
    res <- req(asarDB::make2d(data.table(funtax), dim1 = fl2, dim2 = tl2, numRows = trnum, fun = trfun))


    # if (length(idxM) > 1) {
    #   res <- plotHeatmap(obj, 50, trace = "none", norm = FALSE)
    # } else{
    #   res <- returnAppropriateObj(obj, norm = FALSE, log = TRUE)
    # }
    # res[is.na(res)] <- 0
    numrow1$plot1 <- dim(res)[1]



    if (dim(res)[1] > 1 & dim(res)[2] > 1) {

      downHeat1$is <- TRUE
      plotRes$plot1<-res;
      print(12345)
      res1<<-res
      d3heatmap(
        res,
        dendrogram = chooseDends(res),
        xaxis_height = 220,
        yaxis_width = 270,
        yaxis_font_size = "10px",
        xaxis_font_size = "10px",
        scalecolors = colPal
      )
    } else{
      downHeat1$is <- FALSE
      plotRes$plot1<-NULL;
      cat(555)
      showModal(
        modalDialog(
          title = titleForDimErrorPopup,
          textForDimErrorPopup,
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })


  output$plot11 <- renderPlot({
    cat('output$plot1\n')

    tl1 <- tl1()  # taxonomy selection level


    tl2 <- tl2()  # taxonomy aggregation level
    tn  <- tn()   # taxonomy selection value
    fl1 <- fl1()  # function selection level
    fl2 <- fl2()  # function aggregation level
    fn  <- fn()   # function selection value
    mg1 <- mg1()  # selected sample value

    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation

    cat('plot1-0.1:',trf(),trnum,input$trfun,mg1,'\n')

    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat1$is <- FALSE
      return()
    }
    #  cat(match(fl1,dimNames('tax')),match(fl2,dimNames('tax')))
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))
    cat(paste0('|',c(tl1, tl2, fl1, fl2, mg1),'|'),'\n')#,paste0('|',names(funtaxall),'|'),'\n')
    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
    #    funtax <- funtaxall[, ..keepcols]
    #    funtax <- dplyr::select(funtaxall,c(tl1, tl2, fl1, fl2, mg1))
    #    cat(names(funtax),'\n',tl1, tn, fl1, fn, tl2, fl2,'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2,s=mg1))



    # Apply the selected normalization method to each matrix in the list

    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
      # normalization_factors <- funtax
      # print(normalization_factors)
      #print(333)
    }
    plot_density(as.data.frame(funtax)[,-c(1,2)],norm_method_selected)


  })

  output$text11 <- renderText({


    tl1 <- tl1()  # taxonomy selection level


    tl2 <- tl2()  # taxonomy aggregation level
    tn  <- tn()   # taxonomy selection value
    fl1 <- fl1()  # function selection level
    fl2 <- fl2()  # function aggregation level
    fn  <- fn()   # function selection value
    mg1 <- mg1()  # selected sample value

    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation

    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat1$is <- FALSE
      return()
    }
    #  cat(match(fl1,dimNames('tax')),match(fl2,dimNames('tax')))
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))

    keepcols <- which(names(funtaxall) %in% c(tl1, tl2, fl1, fl2, mg1))
    #    funtax <- funtaxall[, ..keepcols]
    #    funtax <- dplyr::select(funtaxall,c(tl1, tl2, fl1, fl2, mg1))
    #    cat(names(funtax),'\n',tl1, tn, fl1, fn, tl2, fl2,'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2, f2 = fl2,s=mg1))



    # Apply the selected normalization method to each matrix in the list

    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
      # normalization_factors <- funtax
      # print(normalization_factors)
      #print(333)
    }

    paste(names( unlist(summary_stats (as.data.frame(funtax)[,-c(1,2)]))), unlist(summary_stats (as.data.frame(funtax)[,-c(1,2)])),sep=":",collapse = " ")


  })
  ############### Plot1 pdf #########

  plotInput1 <- function(pdf){
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg1 <- mg1()
    colPal <- colPal()
    if(downHeat1$is){
      res<-plotRes$plot1
      main<-makePlotTitle(plot1Title,tl1=tl1,tl2=tl2,tn=tn,fl1=fl1,fl2=fl2,fn=fn,s=mg1,s2=colName())
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          main = 'F vs. T heatmap',
          sepcolor = "black",
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          cex.main=1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if(pdf){
        op = par(mar = c(0, 0, 0, 0),family="mono")
        width<-1100
        height<-1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main<-gsub('\t','    ',main)
        x<-100+strwidth(main,units = 'user',cex = 1)/2
        y<-height-100-strheight(main,units = 'user',cex = 1)/2
        text(x,y,labels = main,adj=0)
        par(op)

      }
    }else{
      return()
    }
  }
  output$downLink1 <- renderUI({
    if(downHeat1$is==TRUE){
      downloadButton(outputId = "down1", label = "Download the heatmap")
    }
  })
  output$down1 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png"){
        png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
        plotInput1(pdf=FALSE)
      } else{
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput1(pdf=TRUE)
      }
      dev.off()
    })

  ############### Heatmap Plot2 #########

  output$plot2 <- renderD3heatmap({
    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation
    cat('plot2-0.1:',trf(),trnum,input$trfun,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat2$is <- FALSE
    }
    classes<-as.character(mdt[match(mg2,sampleNames()),colName()])
    # cat('plot2-1:',tl1,tn,fl1,fl2,fn,mg2,'\n',colName(),sampleNames(),'\n\t mdt',
    #    as.character(mdt[match(mg2,sampleNames()),colName()]),'\n')
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
    # funtax <- funtaxall[,..keepcols]
    # cat(dim(funtax),'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, f2 = fl2,s=mg2,aggSamples=classes))

    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }

    res <- req(asarDB::make2d(data.table(funtax), dim1 = fl2, numRows = trnum, fun = trfun))

    #    funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,f2 = fl2)
    cat(dim(funtax),'\n')
    # obj <- funtax[,-c(1)]
    # cat(dim(obj),colnames(obj),colnames(funtax),'\n')
    # rownames(obj)<-funtax[,fl2]
    # #colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    # cat(colnames(mdt),match(colName(),colnames(mdt)),'\n')
    # cat(dim(obj),colnames(obj),colName(),'\n')
    # dk6 <- data.frame(Means=rowMeans(obj))
    # obj <- obj[which(dk6$Means!=0),]
    # obj <- as.matrix(obj)
    # if(dim(obj)[1]>1){
    #   res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
    # }else{
    #   res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
    # }
    # res[is.na(res)] <- 0
    numrow1$plot2 <- dim(res)[1]
    #main<-makePlot2Title(tl1,tn,fl1,fl2,fn,mg2)
    if(dim(res)[1]>1 & dim(res)[2]>1){
      downHeat2$is <- TRUE
      plotRes$plot2<-res;
      res2<<-res
      d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
    }else{
      downHeat2$is <- FALSE
      showModal(modalDialog(
        title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
      ))
    }
  })

  output$plot21 <- renderPlot({

    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation
    cat('plot2-0.1:',trf(),trnum,input$trfun,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat2$is <- FALSE
    }
    classes<-as.character(mdt[match(mg2,sampleNames()),colName()])
    # cat('plot2-1:',tl1,tn,fl1,fl2,fn,mg2,'\n',colName(),sampleNames(),'\n\t mdt',
    #    as.character(mdt[match(mg2,sampleNames()),colName()]),'\n')
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
    # funtax <- funtaxall[,..keepcols]
    # cat(dim(funtax),'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, f2 = fl2,s=mg2,aggSamples=classes))
    funtax <<-funtax
    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }
    funtax=as.data.frame(funtax)
    funtax=funtax[,-1]
    plot_density(funtax,norm_method_selected)

  })

  output$text21 <- renderText({


    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-trnum()     # number of rows for heatmap truncation
    cat('plot2-0.1:',trf(),trnum,input$trfun,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat2$is <- FALSE
    }
    classes<-as.character(mdt[match(mg2,sampleNames()),colName()])
    # cat('plot2-1:',tl1,tn,fl1,fl2,fn,mg2,'\n',colName(),sampleNames(),'\n\t mdt',
    #    as.character(mdt[match(mg2,sampleNames()),colName()]),'\n')
    validate(
      need(fl1=='toplevel'|(match(fl1,dimNames('fun'))>=
                              match(fl2,dimNames('fun'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, fl1, fl2, mg2))
    # funtax <- funtaxall[,..keepcols]
    # cat(dim(funtax),'\n')
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, f2 = fl2,s=mg2,aggSamples=classes))

    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }

    funtax=as.data.frame(funtax)
    funtax=funtax[,-1]
    paste(names( unlist(summary_stats (funtax))), unlist(summary_stats (funtax)),sep=":",collapse = " ")


  })


  ############### Plot2 #########

  plotInput2 <- function(pdf){
    tl1 <- tl1()
    tn  <- tn()
    fl1 <- fl1()
    fl2 <- fl2()
    fn  <- fn()
    mg2 <- mgall()
    colPal <- colPal()
    if(downHeat2$is){
      res<-plotRes$plot2
      main<-makePlotTitle(plot2Title,tl1=tl1,tn=tn,fl1=fl1,fl2=fl2,fn=fn,s=mg2,s2=colName())
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          sepcolor = "black",
          main=plot2Title,
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          cex.main=1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if(pdf){
        op = par(mar = c(0, 0, 0, 0),family="mono")
        width<-1100
        height<-1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main<-gsub('\t','    ',main)
        x<-100+strwidth(main,units = 'user',cex = 1)/2
        y<-height-100-strheight(main,units = 'user',cex = 1)/2
        text(x,y,labels = main,adj=0)
        par(op)

      }
    }else{
      return()
    }
  }
  output$downLink2 <- renderUI({
    if(downHeat2$is==TRUE){
      downloadButton(outputId = "down2", label = "Download the heatmap")
    }
  })
  output$down2 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png"){
        png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
        plotInput2(FALSE)
      } else{
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput2(TRUE)
      }
      dev.off()
    })


  ############### Heatmap Plot3 #########

  output$plot3 <- renderD3heatmap({
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-as.integer(trnum())     # number of rows for heatmap truncation
    cat('plot3-0.1:',trf(),trnum,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat3$is <- FALSE
    }
    classes<-as.character(mdt[match(mg3,sampleNames()),colName()])
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
    # funtax <- funtaxall[,..keepcols]
    # funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2,s=mg3,aggSamples=classes))

    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }

    res <- req(asarDB::make2d(data.table(funtax), dim1=tl2,numRows = trnum,fun = trfun))

    # obj <- funtax[,-c(1)]
    # rownames(obj)<-funtax[,tl2]
    # #colnames(obj)<-as.character(mdt[c(colnames(obj)), colName()])
    # dk6 <- data.frame(Means=rowMeans(obj))
    # obj <- obj[which(dk6$Means!=0),]
    # obj <- as.matrix(obj)
    # main<-makePlot3Title(tl1,tl2,tn,fl1,fn,mg3)
    # if(dim(obj)[1]>1){
    #   res<-plotHeatmap(obj,50,trace = "none",norm=FALSE)
    # }else{
    #   res<-returnAppropriateObj(obj,norm = FALSE,log = TRUE)
    # }
    # res[is.na(res)] <- 0
    numrow1$plot3 <- dim(res)[1]
    if(dim(res)[1]>1 & dim(res)[2]>1){
      downHeat3$is <- TRUE
      plotRes$plot3<-res;
      res3<<-res
      d3heatmap(res,dendrogram = chooseDends(res), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
    }else{
      downHeat3$is <- FALSE
      showModal(modalDialog(
        title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
      ))
    }
  })



  output$plot31 <- renderPlot({

    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-as.integer(trnum())     # number of rows for heatmap truncation
    cat('plot3-0.1:',trf(),trnum,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat3$is <- FALSE
    }
    classes<-as.character(mdt[match(mg3,sampleNames()),colName()])
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
    # funtax <- funtaxall[,..keepcols]
    # funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2,s=mg3,aggSamples=classes))

    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }
    funtax<<-funtax
    funtax=as.data.frame(funtax)
    funtax=funtax[,-1]
    plot_density(funtax,norm_method_selected)

  })

  output$text31 <- renderText({

    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    trfun   <-getFunction(trf())    # function for heatmap truncation
    trnum   <-as.integer(trnum())     # number of rows for heatmap truncation
    cat('plot3-0.1:',trf(),trnum,'\n')
    colPal <- colPal()
    norm_method_selected <- norm_method()  # Read the selected normalization method
    if(is.null(fn)|is.null(tn)){
      downHeat3$is <- FALSE
    }
    classes<-as.character(mdt[match(mg3,sampleNames()),colName()])
    validate(
      need(tl1=='toplevel'|(match(tl1,dimNames('tax'))>=
                              match(tl2,dimNames('tax'))),
           'The aggregation level should be lower than selection one'))
    # keepcols<-which(names(funtaxall)%in%c(tl1, tl2, fl1, mg3))
    # funtax <- funtaxall[,..keepcols]
    # funtax <- Intfuntax(funtax,tl1,tn,fl1,fn,t2 = tl2)
    funtax <- req(Intfuntax(funtaxall, tl1, tn, fl1, fn, t2 = tl2,s=mg3,aggSamples=classes))

    # Apply the selected normalization method to the data
    if (norm_method_selected != "none") {

      funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
    }
    funtax=as.data.frame(funtax)
    funtax=funtax[,-1]
    paste(names( unlist(summary_stats (funtax))), unlist(summary_stats (funtax)),sep=":",collapse = " ")


  })

  ############### Plot3 pdf #########

  plotInput3 <- function(pdf){
    tl1 <- tl1()
    tl2 <- tl2()
    tn  <- tn()
    fl1 <- fl1()
    fn  <- fn()
    mg3 <- mgall()
    colPal <- colPal()
    if(downHeat2$is){
      res<-plotRes$plot3
      main<-makePlotTitle(plot3Title,tl1=tl1,tl2=tl2,tn=tn,fl1=fl1,fn=fn,s=mg3,s2=colName())
      x <-
        heatmap.2(
          res,
          dendrogram = chooseDends(res),
          col = colPal,
          main=plot3Title,
          sepcolor = "black",
          sepwidth = c(0.05, 0.05),
          key = TRUE,
          keysize = 0.75,
          key.par = list(cex = 0.7),
          symkey = FALSE,
          density.info = "none",
          cexRow = 1,
          cexCol = 1,
          margins = c(20, 30),
          trace = "none",
          srtCol = 50
        )
      if(pdf){
        op = par(mar = c(0, 0, 0, 0),family="mono")
        width<-1100
        height<-1800
        plot(
          c(0, width),
          c(0, height),
          type = "n",
          xlab = "",
          ylab = "",
          xaxs = "i",
          yaxs = "i"
        )
        main<-gsub('\t','    ',main)
        x<-100+strwidth(main,units = 'user',cex = 1)/2
        y<-height-100-strheight(main,units = 'user',cex = 1)/2
        text(x,y,labels = main,adj=0)
        par(op)

      }
    }
  }
  output$downLink3 <- renderUI({
    if(downHeat3$is==TRUE){
      downloadButton(outputId = "down3", label = "Download the heatmap")
    }
  })
  output$down3 <- downloadHandler(
    filename =  function() {
      paste(input$filename, input$var3, sep=".")
    },
    # content is a function with argument file. content writes the plot to the device
    content = function(file) {
      if(input$var3 == "png"){
        png(file, width = 3000, height = 2300, pointsize = 35) # open the png device
        plotInput3(FALSE)
      } else{
        pdf(file, width = 15, height = 15) # open the pdf device
        plotInput3(TRUE)
      }
      dev.off()
    })
  ############### KEGG heatmap #########

  observeEvent(input$path, {
    output$PathwayID <- renderUI({
      tn <- tn()
      if(!is.null(tn)){
        tl1 <- tl1()
        mgall <- mgall()
        ko_sd <- ko_sd()
        trfun   <-getFunction(trf())    # function for heatmap truncation
        trnum   <-trnum()     # number of rows for heatmap truncation
        cat('plot4-0.1:',trf(),trnum,'\n')
        keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
        funtax <- funtaxall[,..keepcols]
        if(tl1=='toplevel'){
          keepcols<-c('toplevel',keepcols)
          funtax$toplevel<-factor('toplevel')
          sp.li<-'toplevel'
          #save(funtax,file='funtax.tmp.Rdata')
        }
        #        funtax <- Intfuntax(funtax,tl1,tn,'toplevel',NULL)
        #        names(funtax)[names(funtax) == tl1] <- 'usp'
        funtax <- Intfuntax(funtaxall,t1=tl1,
                            tn=tn,
                            f1='toplevel',
                            fn='toplevel',
                            t2='md5',
                            f2='ufun',s=mgall)
        #        names(funtax)[names(funtax) == tl1] <- 'usp'
        pathandnames <- as.matrix(getPathwayList(funtax, sp.li =  tn, mgm =  mgall, ko_sd = ko_sd))
        selectInput(inputId = "PathwayID", label = "Input Pathway ID", choices = setNames(as.vector(pathandnames[, "ko"]), pathandnames[,"name"]))
      }})})

  observeEvent(input$goButton, {
    output$plot4 <- renderD3heatmap({
      tl1 <- tl1()
      tn  <- tn()
      mgall <-mgall()
      ko_sd <- ko_sd()
      colPal <- colPal()
      if(is.null(tn)){
        downHeat4$is <- FALSE
        return()
      }
      cat(tl1,tn,mgall,ko_sd,'\n')
      keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
      cat(keepcols,'\n')
      funtax <- funtaxall[,..keepcols]
      names(funtax)[names(funtax) == tl1] <- 'usp'
      cat('plot4:',dim(funtax),'\n',names(funtax))
      obj<-pathwayHeatmap(funtax, tn, mgall, ko_sd)
      # Apply normalization step
      norm_method_selected <- norm_method()  # Read the selected normalization method


      mat3 <- plotHeatmap(obj, 100, normMethod = norm_method_selected, log = TRUE, fun = sd)

      mat3[is.na(mat3)] <- 0
      numrow1$plot4 <- dim(mat3)[1]
      #main<-makePlot4Title(tl1,tn,ko_sd,mgall)
      if(dim(mat3)[1]>1 & dim(mat3)[2]>1){
        downHeat4$is <- TRUE
        d3heatmap(mat3,dendrogram = chooseDends(mat3), xaxis_height = 220, yaxis_width = 270, yaxis_font_size = "10px", xaxis_font_size = "10px", scalecolors = colPal)
      }else{
        downHeat4$is <- FALSE
        showModal(modalDialog(
          title = titleForDimErrorPopup, textForDimErrorPopup, easyClose = TRUE, footer = NULL
        ))
      }
    })})

  ############### KEGG maps #########

  output$Pathway <- renderImage({
    sp.li<- tn()
    tl1 <- tl1()
    pathwi<- pathw()
    mgall <-mgall()
    if(!is.null(sp.li)){
      classes<-as.character(mdt[match(mgall,sampleNames()),colName()])
      # keepcols<-which(names(funtaxall)%in%c(tl1,"ufun","md5", mgall))
      # funtax <- funtaxall[,..keepcols]
      # if(tl1=='toplevel'){
      #   keepcols<-c('toplevel',keepcols)
      #   funtax$toplevel<-factor('toplevel')
      #   sp.li<-'toplevel'
      #   #save(funtax,file='funtax.tmp.Rdata')
      # }
      #funtax <- Intfuntax(funtaxall,t1=tl1,tn=sp.li,f1='toplevel',fn=NULL,f2='ufun',s=mgall,aggSamples=classes)
      funtax <- Intfuntax(funtaxall,t1=tl1,
                          tn=sp.li,
                          f1='toplevel',
                          fn='toplevel',
                          t2='md5',
                          f2='ufun',s=mgall)#,
      #                    aggSamples=classes)
      # Appply normalization
      norm_method_selected <- norm_method()  # Read the selected normalization method
      if (norm_method_selected != "none") {
        funtax <- returnAppropriateObj(funtax, norm = TRUE, log = TRUE, method = norm_method_selected)
      }

      pathImage(funtax, sp.li, mgall, pathwi, kostat, classes)
    }
    cat('src =', paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.leg.png"),'\n')
    list(src = paste0(getwd(),"/","ko", pathwi, ".", sp.li, ".ko.multi.png.multi.png"),
         contentType = 'png',
         alt = "Press GO to select Pathway!")
  }, deleteFile = FALSE)

  ############### Render UIs #########

  output$Mgall <- renderUI({
    selectInput(inputId = "mgall",
                label = metagenomeone,
                choices = sampleNames(),
                selected = c(sampleNames()[metagenome1selected]),
                selectize = TRUE, multiple = TRUE)
  })
  output$mg1 <- renderUI({
    selectInput(inputId = "mg1", label = metagenometwo,
                choices = sampleNames(),
                selected = c(sampleNames()[metagenome2selected]),
                selectize = FALSE)
  })
  output$dynamic1 <- renderUI({
    d3heatmapOutput("plot1", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })
  output$dynamic11 <- renderUI({
    plotOutput("plot11", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })
  output$dynamic111 <- renderUI({
    textOutput("text11")
  })

  output$dynamic2 <- renderUI({
    d3heatmapOutput("plot2", height = paste0(numrow1$plot2*input$pix2+220, "px"))
  })

  output$dynamic21 <- renderUI({
    plotOutput("plot21", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })
  output$dynamic211 <- renderUI({
    textOutput("text21")
  })

  output$dynamic3 <- renderUI({
    d3heatmapOutput("plot3", height = paste0(numrow1$plot3*input$pix3+220, "px"))
  })

  output$dynamic31 <- renderUI({
    plotOutput("plot31", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })
  output$dynamic311 <- renderUI({
    textOutput("text31")
  })

  output$dynamic4 <- renderUI({
    d3heatmapOutput("plot4", height = paste0(numrow1$plot4*input$pix4+220, "px"))
  })

  output$dynamic41 <- renderUI({
    plotOutput("plot41", height = paste0(numrow1$plot1*input$pix1+220, "px"))
  })
  output$dynamic411 <- renderUI({
    textOutput("text41")
  })
  output$taxNames <- renderUI({x <- input$tl1
  if(x!="toplevel"){
    selectInput(inputId = "tn", label = taxthree, multiple=TRUE, choices = as.vector(unique(funtaxall[,get(x)])), selected = taxnames$tn)
  }})
  output$funNames <- renderUI({y <- input$fl1
  if(y!="toplevel"){
    selectInput(inputId = "fn", label = functhree, multiple=TRUE, choices = as.vector(unique(funtaxall[,get(y)])), selected = taxnames$fn)
  }})

  output$colNameO <- renderUI({
    cat('render selector',dim(mdt),'\n')
    selectInput("colName", ColNameSelectorTitle, choices = colnames(rv$data), selected = "colName", selectize = FALSE)
  })

  observeEvent(input$colName,{
    cat('change colName',input$colName,'\n')
    cat(getSampleNames(funtaxall), '\n',rownames(mdt),'\n')
    sampleNames(setNames(getSampleNames(funtaxall), mdt[,input$colName]))
  })

  ############### Palette definitions #########
  output$paletteOutput <- renderPlot({
    display.brewer.pal(8,input$colorPalette)
  }, height = 200, width = 500)
  observeEvent(input$showAllCols,{
    showModal(modalDialog(
      renderPlot({display.brewer.all()}, height = 700, width = 500),title = "All color palettes", easyClose = TRUE, footer = NULL, size = "l"))
  })
  ############### Metadata definitions #########
  output$ui_newcolname <- renderUI({
    textInput("newcolumnname", "Name a new column", sprintf("newcol%s", 1+ncol(values[["DF"]])))
  })
  rv <- reactiveValues(data = mdt)

  output$metadata = renderDT({

    rv$data
  })

  observeEvent(input$loadRdata, {
    inFile <- input$InFile
    cat(class(inFile),inFile$datapath,'\n')
    if(file.exists(inFile$datapath)){
      listA <- loadRdata(inFile$datapath)
      cat(length(listA),'\n')
      dfP<-dim(funtaxall)
      cat(names(funtaxall),'\n')
      funtaxall <<- listA$funtaxall
      cat(dfP,dim(funtaxall),'\n')
      cat(names(funtaxall),'\n')
      mdt <<- listA$mdt
      #colNameS()<-colnames(mdt)
      rv$data<-mdt
      #DF <<- data()
      #input$hot<-mdt
      #values[["hot"]] <- mdt
      d.kres <<- listA$d.kres
      kegg <<- listA$kegg
      ko.path.name <<- listA$ko.path.name
      kostat <<- kostat(funtaxall, d.kres)
      showModal(
        modalDialog(
          title = "Congratulations!",
          "Other Rdata was successfully uploaded",
          easyClose = TRUE,
          footer = NULL
        )
      )
    }
  })

  # pressed <- reactiveValues(a=TRUE)
  #
  # observeEvent(input$addcolumn, {
  #   pressed$a <- FALSE
  #   DF <- isolate(values[["DF"]])
  #   values[["previous"]] <- DF
  #   newcolumn <- eval(parse(text=sprintf('%s(nrow(DF))', isolate(input$newcolumntype))))
  #   values[["DF"]] <- setNames(cbind(DF, newcolumn, stringsAsFactors=FALSE), c(names(DF), isolate(input$newcolumnname)))
  # })
  #
  # observeEvent(input$addcolumn2, {
  #   pressed$a <- TRUE
  # })
  #
  # values = reactiveValues()
  #
  #
  # data = reactive({
  #   if (is.null(input$hot)) {
  #     hot = mdt
  #   } else {
  #     hot = hot_to_r(input$hot)
  #   }
  #
  #   # this would be used as a function input
  #   values[["hot"]] = hot
  #   hot
  # })
  #
  # output$hot <- renderRHandsontable({
  #   if (pressed$a==FALSE){
  #     DF <- values[["DF"]]
  #     if (!is.null(DF)){
  #       rhandsontable(DF, useTypes = FALSE, stretchH = "all")
  #     }
  #   } else {
  #     DF = data()
  #     if (!is.null(DF))
  #       rhandsontable(DF, useTypes = TRUE, stretchH = "all")
  #   }
  # })

})
