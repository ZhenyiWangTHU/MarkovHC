#' Markov Hierarchical Clustering Algorithm
#'
#' Function \code{MarkovHC} The main function of MarkovHC, which applies Markov
#' Hierarchical Clustering Algorithm to the given data (a Seurat object or a
#' matrix/dataframe).
#' @param MarkovHC_input A Seurat object (recommend) or a matrix/dataframe with row names and column names.
#' Each row should be a feature and each column should be a sample.
#' @param SNNslot A character vector indicates the name of the SNN graph slot in the seurat object. If the MarkovHC_input
#' is a matrix, this parameter will be ignored.
#' @param KNNslot A character vector indicates the name of the KNN graph slot in the seurat object. If the MarkovHC_input
#' is a matrix, this parameter will be ignored.
#' @param KNN An interger indicates the number of neighbors in building the
#' KNN graph. Default is 20. If the MarkovHC_input is a Seurat object, the KNN
#' parameter will be set as that in the Seurat object.
#' @param dobasecluster A Bloolean parameter indicates wether to do Louvain clustering
#' on the first level. Default is TRUE.
#' @param cutpoint A numeric value in [0,1] indicates the threshold of
#' the quantile of pseudo energy, default is 0.001.
#' @param verbose A Bloolean parameter indicates wether to print details when run
#'  the program.
#' @details The data given by \code{MarkovHC_input} will be clustered by using
#' the Markov Hierarchical Clustering Algorithm, which generates a hierarchical
#' structure based on the metastability of exponentially perturbed Markov chain.
#' More details of the algorithm please refer to the article.
#' @return This function returns a list including three components, hierarchicalStructure, inputParameters,
#' and midResults, which are the hierarchical clustering results, the input parameters, and some intermediate results.
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export
MarkovHC = function(MarkovHC_input=NULL,
                    SNNslot=NULL,
                    KNNslot=NULL,
                    KNN=20,
                    dobasecluster=TRUE,
                    cutpoint=0.001,
                    verbose = TRUE){

  ##check the input parameters------------------------------------------
  #MarkovHC_input
  if((!is.matrix(MarkovHC_input))&(!(class(MarkovHC_input)=='Seurat'))){
    print("The type of 'MarkovHC_input' must be matirx or seurat object!")
    return(NULL)
  }

  if(is.matrix(MarkovHC_input)){
    print("The input is a matrix.")
    origin_matrix <- t(MarkovHC_input)
  }else{
    print("The input is a Seurat object.")
  }

  #KNN
  if(!is.numeric(KNN)){
    print("The type of 'KNN' must be numeric!")
    return(NULL)
  }

  #cutpoint
  if(!is.numeric(cutpoint)){
    print("The type of 'cutpoint' must be numeric!")
    return(NULL)
  }

  set.seed(1)
  #Do parallel
  ncore<-detectCores()
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)

  ##calculate SNN and build a KNN graph--------------
  if(is.matrix(MarkovHC_input)){
    symmetric_KNN_graph <- FindNeighbors(object = origin_matrix,
                                         k.param = KNN,
                                         compute.SNN = TRUE,
                                         prune.SNN = 0)
    symmetric_KNN_graph <- Seurat::as.sparse(symmetric_KNN_graph$snn)*Seurat::as.sparse(symmetric_KNN_graph$nn)*Matrix::t(Seurat::as.sparse(symmetric_KNN_graph$nn))
  }else{
    symmetric_KNN_graph <- Seurat::as.sparse(MarkovHC_input@graphs[[SNNslot]])*Seurat::as.sparse(MarkovHC_input@graphs[[KNNslot]])*Matrix::t(Seurat::as.sparse(MarkovHC_input@graphs[[KNNslot]]))
  }

  gc(verbose=FALSE)

  ##other options are calculating the degree or PageRank of each node in the graph, we take the degree as the density of the node here.
  #use 'degree' in igraph to calculate the degree of vertexes in the network.

  #use graph_from_adjacency_matrix for convenience
  rownames(symmetric_KNN_graph) <- as.character(1:nrow(symmetric_KNN_graph))
  colnames(symmetric_KNN_graph) <- as.character(1:ncol(symmetric_KNN_graph))
  symmetric_KNN_graph_object <- graph_from_adjacency_matrix(adjmatrix = symmetric_KNN_graph,
                                                            mode = 'undirected',
                                                            weighted = TRUE,
                                                            diag = TRUE)

  centrality_scores <- igraph::degree(symmetric_KNN_graph_object, v = V(symmetric_KNN_graph_object),
                                      mode = "total",
                                      loops = TRUE, normalized = FALSE)

  ##do preclustering----------------------------------------------------
  if(dobasecluster==TRUE){
    #do clustering on the first level
    cluster_louvain_object <- cluster_louvain(graph = symmetric_KNN_graph_object
                                              # weights = E(symmetric_KNN_graph_object)$weight
    )
    hresult_cut <- cluster_louvain_object$memberships[1,]
  }else{
    #do not do clustering on the first level
    hresult_cut <- 1:nrow(symmetric_KNN_graph)
  }

  unique_clusters <- unique(hresult_cut)

  if(dobasecluster==TRUE){
    # merge a cluster as a single point
    # the similarity is the maximum similarity among the cluster group with the other cluster group
    # merge rows
    symmetric_KNN_graph_merged <- as(Matrix::Matrix(data=0, nrow=length(unique_clusters), ncol=ncol(symmetric_KNN_graph), sparse = TRUE), 'dgCMatrix')
    for(clusterindex in unique_clusters){
      temp_index <- which(hresult_cut==clusterindex)
      if(length(temp_index)==1){
        symmetric_KNN_graph_merged[clusterindex,] <- symmetric_KNN_graph[temp_index, ]
      }else{
        temp_cluster <- symmetric_KNN_graph[temp_index, ]
        symmetric_KNN_graph_merged[clusterindex,] <- as.vector(qlcMatrix::colMax(temp_cluster, which = FALSE, ignore.zero = FALSE))
        #symmetric_KNN_graph_merged[clusterindex,] <- apply(temp_cluster, 2, max)
      }
    }

    #merge columns
    #each elements in symmetric_KNN_graph_cluster is the similarity among cluster groups
    symmetric_KNN_graph_cluster <- matrix(0,nrow(symmetric_KNN_graph_merged),nrow(symmetric_KNN_graph_merged))
    rownames(symmetric_KNN_graph_cluster) <- as.character(1:nrow(symmetric_KNN_graph_cluster))
    colnames(symmetric_KNN_graph_cluster) <- as.character(1:ncol(symmetric_KNN_graph_cluster))
    for(clusterindex in unique_clusters){
      for(clusterindex2 in unique_clusters){
        temp_cluster <- symmetric_KNN_graph_merged[clusterindex, which(hresult_cut==clusterindex2)]
        symmetric_KNN_graph_cluster[clusterindex, clusterindex2] <- max(temp_cluster)
      }
    }

    #maximum Jaccard index
    diag(symmetric_KNN_graph_cluster) <- 1

    #calculate the centrality_scores of clusters
    centrality_scores_cluster <- integer(length = length(unique_clusters))
    for (score_index in unique_clusters) {
      centrality_scores_tpm <- centrality_scores[which(hresult_cut==score_index)]
      centrality_scores_cluster[score_index] <- mean(centrality_scores_tpm)
    }
  }else{
    symmetric_KNN_graph_cluster <- as.matrix(symmetric_KNN_graph)
    rownames(symmetric_KNN_graph_cluster) <- as.character(1:nrow(symmetric_KNN_graph_cluster))
    colnames(symmetric_KNN_graph_cluster) <- as.character(1:ncol(symmetric_KNN_graph_cluster))
    centrality_scores_cluster <- centrality_scores
  }

  ## Main part of MarkovHC algorithm
  ## Calculate the transition probability matrix and the pseudo energy matrix
  ## initialize the density matrix, D.matrix
  D.matrix <- matrix(0, nrow = ncol(symmetric_KNN_graph_cluster), ncol = ncol(symmetric_KNN_graph_cluster))
  # set names of rows and columns in D.matrix
  rownames(D.matrix) <- rownames(symmetric_KNN_graph_cluster)
  colnames(D.matrix) <- rownames(symmetric_KNN_graph_cluster)

  nothing <- foreach(nrow_i=1:nrow(D.matrix), .combine='c') %do% {
    D.matrix[nrow_i,] <- (centrality_scores_cluster[nrow_i]/centrality_scores_cluster)^2
    NULL}

  transitionMatrix_pseudoenergyMatrix <- transition_probability_pseudo_energy_matrix(matrix=symmetric_KNN_graph_cluster,
                                                                                     D.matrix=D.matrix)
  #the transition probability matrix
  transitionMatrix<-transitionMatrix_pseudoenergyMatrix$P

  #the pseudo energy matrix
  C_matrix <- transitionMatrix_pseudoenergyMatrix$C

  C_matrix_temp <- C_matrix
  C_matrix_temp[which(abs(C_matrix_temp)<exp(-30))] <- Inf
  ground_energy <- min(C_matrix_temp)/100
  C_matrix <- C_matrix + ground_energy
  C_matrix[which(is.infinite(C_matrix)==TRUE)] <- 0
  rownames(C_matrix) <- as.character(1:nrow(C_matrix))
  colnames(C_matrix) <- as.character(1:ncol(C_matrix))
  C_matrix_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(C_matrix),
                                                       mode = 'directed',
                                                       weighted = TRUE,
                                                       diag = TRUE)

  # Calculate the shortest distance between each vertex pair in the graph
  if(verbose == TRUE){print('Calculate the shortest distance between each vertex pair in the graph.')}

  C_matrix_graph_shortest_distance <- igraph::distances(graph = C_matrix_graph_object,
                                                        v = 1:nrow(C_matrix),
                                                        to = 1:nrow(C_matrix),
                                                        mode = 'out',
                                                        weights = igraph::E(C_matrix_graph_object)$weight,
                                                        algorithm = "dijkstra")

  ## Build the hierarchical structure-----------------------------------
  P_updated <- transitionMatrix
  MarkovHC_result <- list()
  # Store the result of pre-clustering
  attractors <- list()
  basins <- list()
  graphvertex_attractors <- list()
  graphvertex_basins <- list()
  attractorPoints <- list()
  basinPoints <- list()
  for (i in 1:length(unique_clusters)) {
    attractors <- c(attractors, list(i))
    basins <- c(basins, list(i))
    graphvertex_attractors <- c(graphvertex_attractors, list(i))
    graphvertex_basins <- c(graphvertex_basins, list(i))
    clusterPoints <- which(hresult_cut==i)
    attractorPoints <- c(attractorPoints, list(clusterPoints))
    basinPoints <- c(basinPoints, list(clusterPoints))
  }
  basinNum <- length(attractors)
  level_result <- list(basins=basins,
                       attractors=attractors,
                       graphvertex_attractors=graphvertex_attractors,
                       graphvertex_basins=graphvertex_basins,
                       basinPoints=basinPoints,
                       attractorPoints=attractorPoints,
                       basinNum=basinNum,
                       energyCutpoint=0,
                       C_matrix_updatedmatrix=C_matrix_graph_shortest_distance)
  MarkovHC_result <- append(MarkovHC_result, list(level_result))
  levels_indice <- 1
  if(verbose == TRUE){print(paste('Build the level ',as.character(levels_indice),'...', sep = ''))}

  while (TRUE) {
    levels_indice <- levels_indice + 1

    if(verbose == TRUE){print(paste('Build the level ',as.character(levels_indice),'...', sep = ''))}

    ##find basins and attractors
    RS_vector <- judge_RS(P=P_updated)

    ##constructe the list to store the result of this level
    attractors <- list()
    basins <- list()
    graphvertex_attractors <- list()
    graphvertex_basins <- list()
    attractorPoints <- list()
    basinPoints <- list()

    ##partition the state space
    processed_attractors <- integer(length = length(RS_vector))

    processed_attractors[which(RS_vector<=exp(-10))] <- 1
    attractor_indice <- 1
    P_updated_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(P_updated),
                                                          mode = 'directed',
                                                          weighted = TRUE,
                                                          diag = TRUE)

    while(TRUE){
      if(all(processed_attractors==1)){break}
      if(verbose == TRUE){print(paste('Find attractors in the basin ',as.character(attractor_indice),'.', sep = ''))}

      attractor_temp <- which(processed_attractors==0)[1]

      # subcomponent includes the vertex itself in this graph
      attractor_temp_access <- subcomponent(graph = P_updated_graph_object,
                                            v = attractor_temp,
                                            mode = "out")%>%as.vector()
      processed_attractors[attractor_temp_access] <- 1

      attractors <- c(attractors, list(attractor_temp_access))

      basins_temp_merged <- subcomponent(graph = P_updated_graph_object,
                                         v = attractor_temp_access,
                                         mode = "in")%>%as.vector()

      processed_attractors[basins_temp_merged] <- 1

      basins <- c(basins, list(basins_temp_merged))
      attractor_indice <- attractor_indice+1
    }

    ## correct the lost basins caused by R precision limitation
    if(length(unique(unlist(basins)))!=nrow(P_updated)){
      if(verbose == TRUE){print("Correct the lost basins caused by R precision limitation.")}
      lost_basin <- setdiff(1:nrow(P_updated), unique(unlist(basins)))
      lost_attractor <- integer(length = length(lost_basin))
      for(i in 1:length(lost_basin)){
        Si <- subcomponent(graph = P_updated_graph_object,
                           v = lost_basin[i],
                           mode = "out")%>%as.vector()

        Ti <- subcomponent(graph = P_updated_graph_object,
                           v = lost_basin[i],
                           mode = "in")%>%as.vector()
        if(length(setdiff(Si, Ti))>0){
          #transient
          lost_attractor[i] <- 1
        }else{
          #recurrent
          lost_attractor[i] <- 0
        }
      }

      while(TRUE){
        if(all(lost_attractor==1)){break}
        if(verbose == TRUE){print(paste('Find attractors in the basin ',as.character(attractor_indice),'.', sep = ''))}

        attractor_temp <- lost_basin[which(lost_attractor==0)[1]]

        # subcomponent includes the vertex itself in this graph
        attractor_temp_access <- subcomponent(graph = P_updated_graph_object,
                                              v = attractor_temp,
                                              mode = "out")%>%as.vector()
        lost_attractor[which(lost_basin%in%attractor_temp_access)] <- 1

        attractors <- c(attractors, list(attractor_temp_access))

        basins_temp_merged <- subcomponent(graph = P_updated_graph_object,
                                           v = attractor_temp_access,
                                           mode = "in")%>%as.vector()

        lost_attractor[which(lost_basin%in%basins_temp_merged)] <- 1

        basins <- c(basins, list(basins_temp_merged))
        attractor_indice <- attractor_indice+1
      }
    }

    ## assign the points and graph vertexes to basins and attractors
    basinNum <- length(basins)
    basin_indice <- 1
    for (i in 1:basinNum){
      if(verbose == TRUE){print(paste('Partition the basin ',as.character(basin_indice),'.', sep = ''))}

      # assign attractor points
      attractorPoints_temp <- level_result$attractorPoints[attractors[[i]]]%>%unlist()%>%unique()
      attractorPoints <- c(attractorPoints, list(attractorPoints_temp))
      # assign basin points
      basinPoints_temp <- level_result$basinPoints[basins[[i]]]%>%unlist()%>%unique()
      basinPoints <- c(basinPoints, list(basinPoints_temp))
      # assign graphvertex_attractors
      graphvertex_attractors_temp <- level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique()
      graphvertex_attractors <- c(graphvertex_attractors, list(graphvertex_attractors_temp))
      # assign graphvertex_basins
      graphvertex_basins_temp <- level_result$graphvertex_basins[basins[[i]]]%>%unlist()%>%unique()
      graphvertex_basins <- c(graphvertex_basins, list(graphvertex_basins_temp))
      basin_indice <- basin_indice+1
    }

    # update the pseudo energy matrix
    if(verbose == TRUE){print('Update the pseudo energy matrix.')}
    changed_basins <- which((sapply(basins, length)>1)==TRUE)
    unchanged_basins <- which((sapply(basins, length)==1)==TRUE)

    C_matrix_updated <- matrix(data = 0, nrow = basinNum, ncol = basinNum)
    C_matrix_updated[unchanged_basins, unchanged_basins] <- level_result$C_matrix_updatedmatrix[unlist(basins[unchanged_basins]), unlist(basins[unchanged_basins])]

    for (i in changed_basins) {
      for (j in changed_basins) {
        source_v <- level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique()
        target_v <- level_result$graphvertex_basins[basins[[j]]]%>%unlist()%>%unique()
        C_matrix_updated[i,j] <- C_matrix_graph_shortest_distance[source_v, target_v]%>%min()
      }
    }

    for (i in changed_basins) {
      for (j in unchanged_basins) {
        source_v <- level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique()
        target_v <- level_result$graphvertex_basins[basins[[j]]]%>%unlist()%>%unique()
        C_matrix_updated[i,j] <- C_matrix_graph_shortest_distance[source_v, target_v]%>%min()
        #
        source_v <- level_result$graphvertex_attractors[attractors[[j]]]%>%unlist()%>%unique()
        target_v <- level_result$graphvertex_basins[basins[[i]]]%>%unlist()%>%unique()
        C_matrix_updated[j,i] <- C_matrix_graph_shortest_distance[source_v, target_v]%>%min()
      }
    }

    ## update the transition probability matrix

    C_matrix_updated_indice <- C_matrix_updated
    diag(C_matrix_updated_indice) <- Inf
    if((basinNum>1)&(all(is.infinite(C_matrix_updated_indice))==FALSE)){
      if(verbose == TRUE){print('Update the transition probability matrix.')}
      update_P_result <- update_P(C_matrix_updated=C_matrix_updated, C_cut=cutpoint)
      P_updated <- update_P_result[[1]]
      energyCutpoint <- update_P_result[[2]]
    }

    ## constructe the list to store the result of MarkovHC algorithm

    level_result <- list(basins=basins,
                         attractors=attractors,
                         graphvertex_attractors=graphvertex_attractors,
                         graphvertex_basins=graphvertex_basins,
                         basinPoints=basinPoints,
                         attractorPoints=attractorPoints,
                         basinNum=basinNum,
                         C_matrix_updatedmatrix=C_matrix_updated,
                         P_updated=P_updated,
                         energyCutpoint=energyCutpoint,
                         changed_basins=changed_basins
    )
    MarkovHC_result <- append(MarkovHC_result, list(level_result))

    if((basinNum==1)|(all(is.infinite(C_matrix_updated_indice)))){
      ## Output the results---------------------------------------------
      # The input parameters
      inputParameters <- list(
        KNN=KNN,
        dobasecluster=dobasecluster,
        cutpoint=cutpoint
      )
      # The results among the process
      # The sequence of C_cut
      C_cut_seq <- ground_energy

      for(i in 2:length(MarkovHC_result)){
        C_cut_seq <- c(C_cut_seq, MarkovHC_result[[i]][["energyCutpoint"]])
      }

      midResults <- list(
        symmetric_KNN_graph = symmetric_KNN_graph,
        symmetric_KNN_graph_object = symmetric_KNN_graph_object,
        centrality_scores = centrality_scores,
        centrality_scores_cluster = centrality_scores_cluster,
        symmetric_KNN_graph_cluster = symmetric_KNN_graph_cluster,
        transitionMatrix = transitionMatrix,
        C_matrix = C_matrix,
        C_matrix_graph_object=C_matrix_graph_object,
        C_matrix_graph_shortest_distance = C_matrix_graph_shortest_distance,
        C_cut_seq = C_cut_seq,
        ground_energy = ground_energy
      )
      # The MarkovHC object
      MarkovHC_object <- list(
        hierarchicalStructure = MarkovHC_result,
        inputParameters = inputParameters,
        midResults = midResults
      )
      names(MarkovHC_object$hierarchicalStructure) <- paste(rep('level', length(MarkovHC_result)),
                                                            1:length(MarkovHC_result),
                                                            sep = '')
      stopCluster(cl)
      return(MarkovHC_object)
    }
  }
}
