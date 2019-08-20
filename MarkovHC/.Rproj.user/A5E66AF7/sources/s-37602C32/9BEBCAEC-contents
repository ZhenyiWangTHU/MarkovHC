##Extract the results on several levels and arrange them to be compatible with ggplot
#' @export
sankeyResult =function(MarkovObject=NULL,
                       MarkovLevels=NULL){
result.dataframe <- as.data.frame(1:nrow(MarkovObject[["midResults"]][["symmetric_KNN_graph"]]))
colnames(result.dataframe) <- 'index'
for (i in MarkovLevels) {
  result.dataframe.temp <- data.frame(matrix(nrow = MarkovObject$hierarchicalStructure[[i]]$basinPoints%>%unlist()%>%length(),
                                             ncol = 2))
  colnames(result.dataframe.temp) <- c('index', paste('level', as.character(i), sep = ''))
  dataframeRow <- 1
  for (basin in 1:length(MarkovObject$hierarchicalStructure[[i]]$basinPoints)) {
    basinj <- MarkovObject$hierarchicalStructure[[i]]$basinPoints[[basin]]
    for (point in basinj) {
      result.dataframe.temp[dataframeRow,1] <- point
      result.dataframe.temp[dataframeRow,2] <- basin
      dataframeRow <- dataframeRow + 1
    }
  }
  result.dataframe <- merge(result.dataframe, result.dataframe.temp, by='index')
}
return(result.dataframe)
}
