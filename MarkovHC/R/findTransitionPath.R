#' Find the transition path from basinA to basinB
#'
#' Function \code{findTransitionPath} The function finds the samples on the
#' transition path from basinA to basinB.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param level An integer value indicates the customized level.
#' @param basinA An integer value indicates the customized "from" basin.
#' @param basinB An integer value indicates the cusomized "to" basin.
#' @details The function finds the samples on the transition path from basinA to basinB on
#' the customized level.
#' @return  This function returns a list consists of three components:
#' label is a vector, the indices of 1 in this vector indicate the samples on the transition
#' path.
#'
#' pathVertex is the vector contains indices of base-clusters on the transition path. If the
#' dobasecluster parameter of \code{MarkovHC} was set to FALSE, pathVertex contains the sample
#' indices. If the dobasecluster parameter in \code{MarkovHC} was set to TRUE, pathVertex
#' contains the indices of clusters on level 1.
#'
#' shortestPathLength is pseudo energy cost from basinA to basinB.
#'
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
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
