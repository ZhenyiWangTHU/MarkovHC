##Extract the labels of each sample on customized levels
#' @export
fetchLabels =function(MarkovObject=NULL,
                      MarkovLevels=NULL){
  result.dataframe <- matrix(0,
                             nrow = nrow(MarkovObject[["midResults"]][["symmetric_KNN_graph"]]),
                             ncol = length(MarkovLevels))
  result.dataframe_col <- 1
  for (i in MarkovLevels) {
    for (basin in 1:length(MarkovObject$hierarchicalStructure[[i]]$basinPoints)) {
      basinj <- MarkovObject$hierarchicalStructure[[i]]$basinPoints[[basin]]
      result.dataframe[basinj,result.dataframe_col] <- paste(result.dataframe[basinj,result.dataframe_col],'+',basin, sep='')
    }
    result.dataframe_col <- result.dataframe_col+1
  }
  result.dataframe <- as.data.frame(result.dataframe, stringsAsFactors=FALSE)
  for (i in 1:ncol(result.dataframe)) {
    result.dataframe[,i] <- sub("0\\+", "", result.dataframe[,i])
  }
  colnames(result.dataframe) <- paste('level',MarkovLevels,sep = '')
  return(result.dataframe)
}
