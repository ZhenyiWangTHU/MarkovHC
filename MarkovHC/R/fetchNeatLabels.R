##merge small basins to big basins
#' @export
fetchNeatLabels =function(MarkovObject=NULL,
                          level1=NULL,
                          level2=NULL,
                          labelDataFrame=NULL,
                          dm_matrix=NULL,
                          minimumBasin=NULL){

  #
  labels <- labelDataFrame[,level1]
  labelDataFrame <- apply(labelDataFrame, 2, as.numeric)
  allBasinsOnlevel2 <- unique(labelDataFrame[,level2])
  #merge noise basins to its belonged qualified basins in each basin
  for (basin_i in allBasinsOnlevel2) {
    sublabels <- subset(labelDataFrame[,level1], labelDataFrame[,level2]==basin_i)
    freq_labels <- as.data.frame(table(sublabels))
    noise_basins <- as.numeric(as.character(freq_labels[which(freq_labels[,2]<minimumBasin),1]))
    qualified_basins <- as.numeric(as.character(freq_labels[which(freq_labels[,2]>=minimumBasin),1]))

    #If the basin is small, we regard it as noise and merge it to the closest basin
    C_matrix_updated <- MarkovObject[["hierarchicalStructure"]][[level1]][["C_matrix_updatedmatrix"]]
    basinPoints <- MarkovObject[["hierarchicalStructure"]][[level1]][["basinPoints"]]
    #calculate distances bwt noise basins with qualified basins
    noise_basins_to_qualified_basins <- matrix(Inf,nrow(C_matrix_updated),nrow(C_matrix_updated))
    for (noise_basins_i in noise_basins) {
      for (qualified_basins_i in qualified_basins) {
        noise_basins_to_qualified_basins_temp <- dm_matrix[basinPoints[[noise_basins_i]],
                                                           basinPoints[[qualified_basins_i]]]
        noise_basins_to_qualified_basins[noise_basins_i,qualified_basins_i] <- min(noise_basins_to_qualified_basins_temp)
      }
    }
    #let them cannot be recurrent
    for (noise_basins_i in noise_basins) {
      C_matrix_updated[noise_basins_i,noise_basins_i] <- Inf
    }

    #merge noise basins to its "belonged" qualified basins
    for (noise_basins_i in noise_basins) {
      #this noise_basin is an outlier which cannot reach any other qualified basins
      #use dm_matrix
      if(sum(is.infinite(C_matrix_updated[noise_basins_i,qualified_basins]))==(length(qualified_basins))){
        noise_basins_to_qualified_basins[noise_basins_i,setdiff(1:ncol(noise_basins_to_qualified_basins),qualified_basins)] <- Inf
        cloest_basin_index <- which(noise_basins_to_qualified_basins[noise_basins_i,]==min(noise_basins_to_qualified_basins[noise_basins_i,]))
        #C_matrix_updated[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
        labels[which((labels==noise_basins_i)&(labelDataFrame[,level2]==basin_i))] <- cloest_basin_index[1]
      }else{
        #this noise_basin is not an outlier
        #use C_matrix
        C_matrix_updated[noise_basins_i,setdiff(1:ncol(C_matrix_updated),qualified_basins)] <- Inf
        cloest_basin_index <- which(C_matrix_updated[noise_basins_i,]==min(C_matrix_updated[noise_basins_i,]))
        #C_matrix_updated[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
        labels[which((labels==noise_basins_i)&(labelDataFrame[,level2]==basin_i))] <- cloest_basin_index[1]
      }
    }
  }
  return(labels)

  # labelDataFrame <- apply(labelDataFrame, 2,as.numeric)
  # labels <- labelDataFrame[,level1]
  # freq_labels <- as.data.frame(table(labels))
  # noise_basins <- as.numeric(as.character(freq_labels[which(freq_labels[,2]<minimumBasin),1]))
  # qualified_basins <- as.numeric(as.character(freq_labels[which(freq_labels[,2]>=minimumBasin),1]))
  # print(qualified_basins)
  # #If the basin is small, we regard it as noise and merge it to the closest basin
  # C_matrix_updated <- MarkovObject[["hierarchicalStructure"]][[level1]][["C_matrix_updatedmatrix"]]
  # basinPoints <- MarkovObject[["hierarchicalStructure"]][[level1]][["basinPoints"]]
  # #calculate distances bwt noise basins with qualified basins
  # noise_basins_to_qualified_basins <- matrix(Inf,nrow(C_matrix_updated),nrow(C_matrix_updated))
  # for (noise_basins_i in noise_basins) {
  #   for (qualified_basins_i in qualified_basins) {
  #     noise_basins_to_qualified_basins_temp <- dm_matrix[basinPoints[[noise_basins_i]],
  #                                                        basinPoints[[qualified_basins_i]]]
  #     noise_basins_to_qualified_basins[noise_basins_i,qualified_basins_i] <- min(noise_basins_to_qualified_basins_temp)
  #   }
  # }
  #
  # #row_min <- apply(C_matrix_updated, 1, min)
  # #let them cannot be recurrent
  # for (noise_basins_i in noise_basins) {
  #   C_matrix_updated[noise_basins_i,noise_basins_i] <- Inf
  # }
  # #merge noise basins to its "belonged" qualified basins
  # for (noise_basins_i in noise_basins) {
  #   #this noise_basin is an outlier which cannot reach any other qualified basins
  #   #use dm_matrix
  #   level2qulified_basin <- unique(labelDataFrame[which(labelDataFrame[,level1]==noise_basins_i), level2])
  #
  #   qualified_basins_temp <- intersect(qualified_basins,
  #                                      labelDataFrame[which(labelDataFrame[,level2]==level2qulified_basin),level1])
  #   print(qualified_basins_temp)
  #   if(sum(is.infinite(C_matrix_updated[noise_basins_i,qualified_basins_temp]))==(length(qualified_basins_temp))){
  #     noise_basins_to_qualified_basins[noise_basins_i,setdiff(1:ncol(noise_basins_to_qualified_basins),qualified_basins_temp)] <- Inf
  #     cloest_basin_index <- which(noise_basins_to_qualified_basins[noise_basins_i,]==min(noise_basins_to_qualified_basins[noise_basins_i,]))
  #
  #     print('!')
  #     print(min(noise_basins_to_qualified_basins[noise_basins_i,]))
  #     print(cloest_basin_index)
  #     #C_matrix_updated[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
  #     labels[which(labels==noise_basins_i)] <- cloest_basin_index[1]
  #   }else{
  #     #this noise_basin is not an outlier
  #     #use C_matrix
  #     C_matrix_updated[noise_basins_i,setdiff(1:ncol(C_matrix_updated),qualified_basins_temp)] <- Inf
  #     cloest_basin_index <- which(C_matrix_updated[noise_basins_i,]==min(C_matrix_updated[noise_basins_i,]))
  #     print('!!')
  #     print(min(C_matrix_updated[noise_basins_i,]))
  #     print(cloest_basin_index)
  #     #C_matrix_updated[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
  #     labels[which(labels==noise_basins_i)] <- cloest_basin_index[1]
  #   }
  # }

}
