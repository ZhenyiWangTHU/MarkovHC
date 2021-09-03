#' Select top principle components
#'
#' Function \code{PC_selection} The function selects top the principle components by
#' finding the change point of standard deviations of each principle component.
#' @param SeuratObject A seurat object after RunPCA run.
#' @return  This function plots the standard deviations of the principle components and returns
#' the location(change point)+1.
#' @details The piecewise.linear function in "SiZer" package is used to fit a degree 1 spline
#' with 1 change point and return the location(change point)+1.
#'
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export
#'
PC_selection = function(SeuratObject=NULL){
  model <- SiZer::piecewise.linear(x=1:length(SeuratObject@reductions$pca@stdev), y=SeuratObject@reductions$pca@stdev)
  recommend_level <- ceiling((model$change.point)+1)
  plot(model)
  print(recommend_level)
}
