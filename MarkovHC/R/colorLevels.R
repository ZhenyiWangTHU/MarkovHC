##Customize the colors on each level.
#' @export
colorLevels =function(labels=NULL,
                      colorVector=NULL,
                      basinOrders=NULL){
  color_dataframe <- labels
  color_dataframe[,1] <-  plyr::mapvalues(x = labels[,1],
                                          from = basinOrders,
                                          to = colorVector,
                                          warn_missing = FALSE)
  for (i in 2:ncol(labels)) {
    for (j in unique(labels[,i])) {
      color_dataframe_temp <- color_dataframe[which(labels[,i]==j),c(i-1,i)]
      maxcolor_temp <- as.data.frame(table(color_dataframe_temp[,1]), stringsAsFactors=FALSE)
      maxcolor <- maxcolor_temp[which(maxcolor_temp[,2]==max(maxcolor_temp[,2]))[1],1]
      color_dataframe[which(labels[,i]==j),i] <- maxcolor
    }
  }
  return(color_dataframe)
}

