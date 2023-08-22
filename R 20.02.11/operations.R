#' Aggregation operation for hierarchical dimensions (taxonomy and function).
#'
#' Function takes 3D matrix in the data.table form and aggregate
#'
#' @param matrix -- data.table representing 3D matrix
#' @param aggDim -- name of the dimension. Should be one of 'taxonomy','function'.
#' @param aggLevel -- name of the level at which aggregation is performed.
#'
#' @return 3D data.table aggregated at particular level
#' @export
#'
aggregationH<-function(matrix,aggDim,aggLevel){

}

#' Aggregation operation for non-hierarchical samples dimension.
#'
#' Function requires additional data.frame for metadata to choose column from.
#' Identical values in the selected column will be aggregated.
#'
#' @param matrix -- data.table representing 3D matrix
#' @param metadata -- data.frame describing samples
#' @param aggLevel -- name of the column in the metadata to use for aggregation
#'
#' @return -- data.table representing aggregated 3D matrix
#' @export
#'
aggregationS<-function(matrix,metadata,aggLevel){

}

#' Selection operation for hierarchical dimensions (taxonomy and function).
#'
#' Function takes 3D matrix, name of the dimension, level of the hierarchy and
#' list of values to take into account
#'
#' @param matrix -- data.table representing 3D matrix
#' @param selDim -- name of the dimension. Should be one of 'taxonomy','function'.
#' @param selLevel -- name of the level at which selection is performed.
#' @param selValue -- vector of values to keep in the dataset.
#'
#' @return data.table representing subset of the 3D matrix
#' @export
#'
selectionH<-function(matrix,selDim,selLevel,selValue){

}

#' Selection operation for non-hierarchical samples dimension.
#'
#' Function requires additional data.frame for metadata to choose column from.
#' Only samples described by chosen values in the selected column will be kept.
#'
#' @param matrix -- data.table representing 3D matrix
#' @param metadata -- data.frame describing samples
#' @param selLevel -- name of the column in the metadata to use for
#' @param selValue -- vector of values to keep in the dataset.
#'
#' @return
#' @export
#'
selectionS<-function(matrix,metadata,selLevel,selValue){

}