#' Recommend levels by integrating internal measures
#'
#' Function \code{IMI_selection} The function recommends levels
#' by integrating internal measures.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param prune A Bloolean parameter indicates wether to merge small basins into
#' big basins. Default is TRUE.
#' @param weed An integer value defines small basins. When the prune parameter is
#' set to TRUE, the basins include less points than weed will be merged into
#' big basins on each level.
#' @return  This function returns a dataframe containing the ranked levels and internal measures.
#' @details IMI_selection deploys connectivity, silhouette and Dunn from "clValid" package.
#' And aggregate them with energy gaps in MarkovObject to rank the levels.
#'
#' Refer to the manuscript of "clValid" package, the connectivity has a value between zero and
#' infinite and should be minimized, well-clustered observations having silhouette values near 1
#' and poorly clustered observations having silhouette values near 1, and the Dunn Index has a value
#' between zero and infinite, and should be maximized. Meanwhile, the C_cut value of each level in MarkovHC
#' indicates the energy level, so the first order difference of the C_cut values, which we termed energy gaps(C_cut_gap)
#' indicates the energy difference between neighbor levels and this value should be maximized. We rank levels
#' by these four measures and aggrate the ranks using aggregateRanks function in "RobustRankAggreg" package.
#'
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export

IMI_selection = function(MarkovObject=NULL,
                         prune=NULL,
                         weed=NULL){
  C_matrix_graph_shortest_distance_for_level_selection <- MarkovObject$midResults$C_matrix_graph_shortest_distance
  C_matrix_graph_shortest_distance_for_level_selection[which(is.infinite(C_matrix_graph_shortest_distance_for_level_selection)==TRUE)] <- (-1)
  C_matrix_graph_shortest_distance_for_level_selection[which((C_matrix_graph_shortest_distance_for_level_selection<0)==TRUE)] <- max(C_matrix_graph_shortest_distance_for_level_selection)
  C_matrix_graph_shortest_distance_for_level_selection <- (C_matrix_graph_shortest_distance_for_level_selection+t(C_matrix_graph_shortest_distance_for_level_selection))/2
  labels <-  fetchLabels(MarkovObject=MarkovObject,
                         MarkovLevels=1:length(MarkovObject$hierarchicalStructure),
                         prune = FALSE, weed = 10)
  for(i in 2:ncol(labels)){
    max.temp <- max(table(labels[,i]))
    if(max.temp <= weed){
      next
    }else{
      startLv <- i
      break
      }
  }
  labels_temp <-  fetchLabels(MarkovObject=MarkovObject,
                              MarkovLevels=startLv:length(MarkovObject$hierarchicalStructure),
                              prune = prune, weed = weed)
  labels[,startLv:ncol(labels)] <- labels_temp[,1:ncol(labels_temp)]
  labels_unique<-unique(labels)

  labels_unique <- labels_unique[order(as.numeric(labels_unique[,1]), decreasing = FALSE),]

  #connectivity
  connectivity_levels <- numeric(length = ncol(labels_unique))

  for(i in 1:ncol(labels)){
    connectivity_levels[i] <- connectivity(distance = C_matrix_graph_shortest_distance_for_level_selection,
                                           clusters = as.integer(factor(labels_unique[,i])),
                                           Data = NULL,
                                           neighbSize = ceiling(dim(C_matrix_graph_shortest_distance_for_level_selection)[1]/5))
  }
  connectivity_levels <- abs(diff(connectivity_levels))
  connectivity_levels <- c(min(connectivity_levels), connectivity_levels)
  if(length(unique(labels_unique[,ncol(labels_unique)]))==1){
    connectivity_levels[length(connectivity_levels)] <- min(connectivity_levels)
  }

  #silhouette
  silhouette_levels <- numeric(length =ncol(labels))
  silhouette_levels[1] <- (-1)
  for(i in 2:ncol(labels)){
    if(length(unique(labels_unique[,i]))==1){
      silhouette_levels[i] <- (-1)
    }else{
      silhouette_levels[i] <- summary(silhouette(x=as.numeric(factor(labels_unique[,i])),
                                                 dmatrix=C_matrix_graph_shortest_distance_for_level_selection))$avg.width
    }
  }

  #Dunn
  dunn_levels <- numeric(length = ncol(labels))
  dunn_levels[1] <- 0
  for(i in 2:ncol(labels)){
    if(length(unique(labels_unique[,i]))==1){
      silhouette_levels[i] <- 0
    }else{
      dunn_levels[i] <- dunn(distance = C_matrix_graph_shortest_distance_for_level_selection,
                             clusters = as.integer(factor(labels_unique[,i])),
                             Data = NULL)
    }
  }

  #C_cut_gap
  C_cut_gap <- abs(diff(MarkovObject$midResults$C_cut_seq))
  C_cut_gap <- c(0,C_cut_gap)

  # robust rank aggregation
  level_rank <- list(as.character(order(connectivity_levels, decreasing = TRUE)),
                     as.character(order(silhouette_levels, decreasing = TRUE)),
                     as.character(order(dunn_levels, decreasing = TRUE)),
                     as.character(order(C_cut_gap, decreasing = TRUE)))
  level_rank_results <- RobustRankAggreg::aggregateRanks(glist=level_rank)
  level_rank_results[,1] <- as.integer(as.character(level_rank_results[,1]))
  level_rank_results$connectivity <- mapvalues(level_rank_results[,1], from=1:length(connectivity_levels),to=connectivity_levels)
  level_rank_results$silhouette <- mapvalues(level_rank_results[,1], from=1:length(silhouette_levels),to=silhouette_levels)
  level_rank_results$dunn <- mapvalues(level_rank_results[,1], from=1:length(dunn_levels),to=dunn_levels)
  level_rank_results$C_cut_gap <- mapvalues(level_rank_results[,1], from=1:length(C_cut_gap),to=C_cut_gap)
  return(level_rank_results)
}
