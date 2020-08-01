##Find the transition path from basinA to basinB
#' @export
findTransitionPath = function(MarkovObject = NULL,
                              level = NULL,
                              basinA = NULL,
                              basinB = NULL){
  C_matrix_graph_object <- MarkovObject[["midResults"]][["C_matrix_graph_object"]]
  paths <- list()
  for (i in MarkovObject$hierarchicalStructure[[level]]$graphvertex_attractors[[basinA]]) {
    paths_temp <- all_shortest_paths(graph = C_matrix_graph_object,
                                     from = i,
                                     to = MarkovObject$hierarchicalStructure[[level]]$graphvertex_basins[[basinB]],
                                     mode = 'out',
                                     weights = E(C_matrix_graph_object)$weight)
    paths <- c(paths, paths_temp[[1]])
  }
  if(length(paths)==0){return(list())}
  shortestPathLength <- Inf
  pathVertex <- NULL
  for (i in 1:length(paths)) {
    #the minmum continuous number is always exclusive
    pathLength <- distances(C_matrix_graph_object,
                            v=paths[[i]][1],
                            to=paths[[i]][length(paths[[i]])],
                            mode = 'out',
                            weights = E(C_matrix_graph_object)$weight,
                            algorithm = "dijkstra")
    if(pathLength < shortestPathLength){
      shortestPathLength <- pathLength
      pathVertex <- as.vector(paths[[i]])
    }
  }
  label <- integer(dim(MarkovObject[["midResults"]][["symmetric_KNN_graph"]])[1])
  for (i in pathVertex) {
    label[MarkovObject$hierarchicalStructure[[1]]$basinPoints[[i]]] <- 1
  }
  return(list(label,pathVertex,shortestPathLength))
}
