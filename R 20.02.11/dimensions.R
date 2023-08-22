#' Get dimension definition.
#'
#' Function prepare list of data.table columns indices that constitute selected
#' dimension. In the case that matrix was already aggregated, which means some
#' dimensions are missing, only subset of dimensions that exists in the matrix is
#' returned.
#'
#' @param matrix -- data.table representing 3D matrix
#' @param dim -- name of the dimension. Should be one of 'taxonomy','functions','samples'.
#'
#' @return vector of dimension column indices
#' @export
#'
getDim<-function(matrix,dim=c('taxonomy','functions','samples')){
  dim <- match.arg(dim)
  dNames<-switch(dim,
         taxonomy = dimNames(dim),
         functions = dimNames(dim),
         samples = getSampleNames(matrix))
  n<-colnames(matrix)
  ind<-match(dNames,n)
  return(ind[!is.na(ind)])
}

#' Get dimension names for samples in particular matrix.
#'
#' @param matrix -- data.table representing 3D matrix
#'
#' @return character vector of colunm names
#' @export
#'
getSampleNames<-function(matrix){
  t<-dimNames('tax')
  f<-dimNames('fun')
  exclude<-c(t,f,'md5')
  idx<-which(!colnames(matrix)%in%exclude)
  return(colnames(matrix)[idx])
}

#' Return list of predefined names for taxonomy and function dimensions.
#'
#' @param dim name of dimension should be either 'taxonomy' or 'functions'.
#'
#' @return vector of column names.
#' @export
#'
dimNames<-function(dim=c('taxonomy','functions')){
  dim <- match.arg(dim)
  switch(dim,
         taxonomy = c('usp','species','genus','family','order','class','phylum','domain','kingdom'),
         functions = c('ufun','FUN1','FUN2','FUN3','FUN4'))
 }


#' Get level column index.
#'
#' Function finds the index of the column, corresponding to the particular
#' level at particular hierarchical dimension in the 3D matrix. If
#'
#' @param matrix -- data.table representing 3D matrix
#' @param dim -- name of the dimension. Should be one of 'taxonomy','function'.
#' @param level -- name of the level in the dimension.
#'
#' @return integer index of the column of interest, or NA.
#' @export
#'
getLevel<-function(matrix,dim,level){

}