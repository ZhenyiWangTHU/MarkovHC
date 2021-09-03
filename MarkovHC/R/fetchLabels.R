#' Extract the clustering labels
#'
#' Function \code{fetchLabels} The function extracts the clustering results
#' of each sample on customized levels.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param MarkovLevels An integer vector indicates the customized levels.
#' @param prune A Bloolean parameter indicates wether to merge small basins into
#' big basins. Default is TRUE.
#' @param weed An integer value defines small basins. When the prune parameter is
#' set to TRUE, the basins include less points than weed will be merged into
#' big basins on each level.
#' @details This function extracts the clustering results of each sample on customized levels
#' and returns a dataframe.
#' @return This function returns a dataframe, whose columns correspond to customized levels,
#' rows correspond to samples, and entries are clustering labels.
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export
fetchLabels =function(MarkovObject=NULL,
                      MarkovLevels=NULL,
                      prune=TRUE,
                      weed=10){
  result.dataframe <- data.frame(matrix("A",
                                        nrow = length(MarkovObject[["midResults"]][["centrality_scores"]]),
                                        ncol = length(MarkovLevels)))
  result.dataframe = as.data.frame(lapply(result.dataframe, as.character), stringsAsFactors=FALSE)

  result.dataframe_col <- 1
  for (i in MarkovLevels) {
    for (basin in 1:length(MarkovObject$hierarchicalStructure[[i]]$basinPoints)) {
      basinj <- MarkovObject$hierarchicalStructure[[i]]$basinPoints[[basin]]
      result.dataframe[basinj,result.dataframe_col] <- as.character(paste(result.dataframe[basinj,result.dataframe_col],'+',as.character(basin), sep=''))
    }
    result.dataframe_col <- result.dataframe_col+1
  }
  result.dataframe <- as.data.frame(result.dataframe, stringsAsFactors=FALSE)
  for (i in 1:ncol(result.dataframe)) {
    result.dataframe[,i] <- sub("A\\+", "", result.dataframe[,i])
  }
  colnames(result.dataframe) <- paste('lv',MarkovLevels,sep = '')
  #merge the small basin to its cloest largest isolated basin
  if(prune==TRUE){
    for(i in 1:ncol(result.dataframe)){
      C_matrix_i <- MarkovObject$hierarchicalStructure[[paste('level',MarkovLevels[i],sep = '')]]$C_matrix_updatedmatrix
      basin_size <- table(result.dataframe[,i])
      small_basins <- names(basin_size[which(basin_size<=weed)])
      big_basins <- as.numeric(setdiff(as.character(1:nrow(C_matrix_i)), small_basins))
      for(j in small_basins){
        small_basins_temp <- as.numeric(unlist(str_split(string=j, pattern = '\\+')))
        C_matrix_i_temp <- matrix(C_matrix_i[small_basins_temp, big_basins],nrow = length(small_basins_temp), ncol = length(big_basins))
        cloest_index <- as.data.frame(which(C_matrix_i_temp==min(C_matrix_i_temp), arr.ind = TRUE))$col
        big_basins_temp <- big_basins[cloest_index]
        basin_size_temp <- basin_size[which(names(basin_size)%in%big_basins_temp)]
        cloest_big_basin <- names(which.max(basin_size_temp))
        result.dataframe[which(result.dataframe[,i]==j),i] <- cloest_big_basin
      }
    }
  }

  return(result.dataframe)
}
