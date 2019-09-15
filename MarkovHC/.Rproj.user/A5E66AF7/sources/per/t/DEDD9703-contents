#order points according to the merge order of hierarchical structure
#' @export
orderMarkovHC = function(MarkovObject=NULL,
                         MarkovLevels=NULL,
                         orderLevels=NULL){
  result.dataframe <- as.data.frame(1:nrow(MarkovObject[["midResults"]][["symmetric_KNN_graph"]]))
  colnames(result.dataframe) <- 'index'
  for (i in MarkovLevels) {
    result.dataframe.temp <- data.frame(matrix(nrow = nrow(result.dataframe),
                                               ncol = 2))
    colnames(result.dataframe.temp) <- c('index', paste('level', as.character(i), sep = ''))
    result.dataframe.temp$index <- 1:nrow(result.dataframe.temp)
    for (basin in 1:length(MarkovObject$hierarchicalStructure[[i]]$basinPoints)) {
      basinj <- MarkovObject$hierarchicalStructure[[i]]$basinPoints[[basin]]
      result.dataframe.temp[basinj,2] <- basin
    }
    result.dataframe <- merge(result.dataframe, result.dataframe.temp, by='index')
  }
  orderLevels <- orderLevels[order(orderLevels, decreasing = TRUE)]%>%
                 paste('level',., sep='')%>%
                 paste(., collapse = '+')%>%
                 paste('~',., sep = '')
  result.dataframe <- orderBy(as.formula(orderLevels), result.dataframe)
  return(result.dataframe)
}
