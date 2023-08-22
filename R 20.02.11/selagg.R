#' Selection and aggregation function
#'
#' @param data -- data.table of unprocessed abundances
#' @param taxS  -- taxonomy selection level
#' @param taxVal  -- taxonomy selection value
#' @param funS  -- function selection level
#' @param funVal  -- function selection level
#' @param samples  -- samples to work with
#' @param taxAgg  -- taxonomy aggregation level
#' @param funAgg  -- function aggregation level
#' @param aggSamples  -- list of samples metadata. All samples with the
#' same metadata will be aggregated
#' @param aggFunction  -- aggregation function (default sum)
#'
#' taxS and funS could accomodate value 'top', which means that the selection is not required.
#'
#' @return -- data.table filtered and aggregated
#' @export
#'
selagg <- function(data, taxS, taxVal, funS, funVal,samples, taxAgg=NULL, funAgg=NULL,aggSamples=c(),aggFunction=sum){
  if(length(aggSamples)>0&length(aggSamples)!=length(samples)){
    stop('aggSamples should be either empty, or has the same length as samples')
  }

  if(!"data.table" %in% class(data)){
    data<-data.table(data)
  }
  keepcols <- which(names(data) %in% c(taxS, taxAgg, funS, funAgg, samples))
  result2 <- data[, ..keepcols]
  if(length(taxVal)>1){
    grepTax<-paste0('(',paste0(taxVal,collapse = '|'),')')
  }else{
    grepTax<-taxVal
  }

  if(length(funVal)>1){
    grepFun<-paste0('(',paste0(funVal,collapse = '|'),')')
  }else{
    grepFun<-funVal
  }

  cat('selagg-1:','taxS=',taxS, 'taxVal=',taxVal, 'funS=',funS, 'funVal=',funVal, 'taxAgg=',taxAgg, 'funAgg=',funAgg,'\n')
  cat('selagg-2:',dim(result2),'\n',names(result2),'\n')
  cat('selagg-2.5:',grepTax,grepFun,'\n',class(result2),'\n')
  if(taxS!="toplevel"){
    result2<-dplyr::slice_(result2,paste0("grep('",grepTax,"',",taxS,')'))
  }
  if(funS!="toplevel"){
    result2<-dplyr::slice_(result2,paste0("grep('",grepFun,"',",funS,')'))
  }
  if(length(aggSamples)>0){
    sCols<-which(names(result2) %in% c(samples))
    aggS<-unique(aggSamples)
    cat('selagg-2.8:',aggS,'\n')
    result2<-data.table(result2)
    result2[,(aggS):=0]
    sapply(aggS,function(.x)samples[which(aggSamples==.x)])->idx
    #cat('selagg-2.9:',idx,'\n')
    for(i in 1:length(idx)){
      result2[,c(names(idx)[i]):=rowSums(.SD),.SDcols=c(idx[[i]])]
      result2[,c(idx[[i]]):=NULL]
    }
  }
  if(!is.null(taxAgg) & !is.null(funAgg)){
    result2<- plyr::ddply(result2, c(taxAgg,funAgg), numcolwise(aggFunction))
  }else{
    if(!is.null(taxAgg)){
      result2<- plyr::ddply(result2, taxAgg, numcolwise(aggFunction))
    }
    if(!is.null(funAgg)){
      result2<- plyr::ddply(result2, funAgg, numcolwise(aggFunction))
    }
  }
  cat('selagg-3:',dim(result2),'\n',names(result2),'\n')
  if(any(dim(result2)==0)){
    return(NULL)
  }
  idx<-which(!is.na(result2[,1])& !is.na(result2[,2])& result2[,2]!= "")
  result2 <- data.table(result2[idx,])
  cat('selagg-4:',dim(result2),'\n',names(result2),'\n')
  return(result2)
}


make2d <- function(data,dim1,dim2=NULL,numRows=50,sub.na=0,norm=FALSE,log=TRUE,fun=sd){
  f<-data
  sidx<-getSampleNames(data)
  names(f)<-sub(dim1,'dim1',names(f))
  if(!is.null(dim2)){
    names(f)<-sub(dim2,'dim2',names(f))
    f[,val:=rowSums(.SD),.SDcols=sidx]
    fcast<-dcast(data = f,dim1~dim2,fun.aggregate = sum,value.var = 'val')
  # obj <- matrix(nrow = length(unique(f[,2])), ncol = length(unique(f[,1])))
  # colnames(obj) <- unique(f[,1])
  # rownames(obj) <- unique(f[,2])
  # for (x in 1:nrow(funtax)){
  #   obj[as.character(f[x,2]), as.character(f[x,1])]<- f[x,3]
  # }
  # obj[is.na(obj)]<-sub.na
  #
  # rowmean <- rowMeans(obj)
  # colmean <- colMeans(obj)
  # obj <- obj[which(rowmean$Means != 0), which(colmean$Means != 0)]
  }else{
    fcast<-f[,lapply(.SD,sum),by=dim1,.SDcols=sidx]
  }
  obj<-as.matrix(fcast[,-c(1)])
  rownames(obj)<-fcast$dim1
  return(prepareHeatmap(obj,n=numRows,norm=norm,log=log,fun=fun))
}

#' Prepare data for heatmap
#'
#' @param obj original matrix
#' @param n number of rows to keep
#' @param m maximum number of columns to show
#' @param norm should data be normalized
#' @param log should log2 transformation be applied
#' @param fun the function to use for column and row selections
#' @param ...
#'
#' @return heatmap matrix
#' @export
prepareHeatmap<-function(obj,n,m=50,norm=TRUE,log=TRUE,fun=sd,...){
  mat = normAndScale(obj, norm, log)
  cat('prepareHeatmap-1:',dim(mat),class(mat),'\n')
  colsToKeep = which(colSums(mat,na.rm = TRUE) > 0)
  otusToKeep = which(rowSums(mat,na.rm = TRUE) > 0)
  colIndices = colsToKeep
  cat('prepareHeatmap-2:',otusToKeep,'\n',colsToKeep,'\n')
  if(length(otusToKeep)<=1|length(colsToKeep)<=1){
    return(matrix(-1,ncol = 1,nrow = 1))
  }
  otuStats = apply(mat[otusToKeep, colsToKeep], 1, fun)
  colsStats = apply(mat[otusToKeep, colsToKeep], 2, fun)
  colIndices = colsToKeep[order(colsStats, decreasing = TRUE)[1:min(c(length(colsToKeep),m,dim(mat)[2]))]]
  cat('prepareHeatmap-2.1:',c(length(colsToKeep),m,dim(mat)[2]),'\n',min(c(length(colsToKeep),m,dim(mat)[2])),'\n')
  otuIndices = otusToKeep[order(otuStats, decreasing = TRUE)[1:min(c(length(otusToKeep),n,dim(mat)[1]))]]
  cat('prepareHeatmap-2.2:',c(length(otusToKeep),n,dim(mat)[1]),'\n',min(c(length(otusToKeep),n,dim(mat)[1])),'\n')
  mat2 = mat[otuIndices, colIndices]
  if(dim(mat)[1]>length(otuIndices)&dim(mat)[2]>length(colIndices)){
    mes<- 'Only 50 rows and columns with higher dispersion is shown.'
    attr(mat2,'message')<-mes
  }else if(dim(mat)[1]>length(otuIndices)){
    mes<-'Only 50 rows with higher dispersion is shown.'
    attr(mat2,'message')<-mes
  }else if(dim(mat)[2]>length(colIndices)){
    mes<-'Only 50 columns with higher dispersion is shown.'
    attr(mat2,'message')<-mes
  }

  cat('prepareHeatmap-3:',dim(mat2),'\n')
  return(mat2)
}
normAndScale <- function(obj, norm, log) {
  # cat("normAndScale: dim(obj)",dim(obj),'\n',class(obj),'\n',colnames(obj),'\n')
  # cat("normAndScale: any(dim(obj) == 0)",any(dim(obj) == 0),'\n')
  if (!inherits(obj,'matrix') ||
      any(dim(obj) == 0)){
    #return(matrix(-1, ncol = 1, nrow = 1))
    stop('Data should be a matrix with all dimensions longer than 1')
  }
  res <- avearrays(obj)

  if (log) {
    res <- log2(res + 1)
  }
  if (norm) {
    res <- scale(res)
  }
  return(res)
}
