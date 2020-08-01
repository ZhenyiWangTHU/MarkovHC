#' Markov Hierarchical Clustering Algorithm
#'
#' Function \code{MarkovHC} The main function of MarkovHC, which applies Markov
#' Hierarchical Clustering Algorithm to the given data (a dataframe or a matrix).
#' @param origin_matrix A matrix or a dataframe with row names and column names.
#' Each row should be a feature and each column should be a sample.
# #' @param minrt A parameter helping to judge whether a group of samples is a
# #' qualified cluster. If we assume that we have \code{tol_num} samples in
# #' total, then a qualified cluster should be of size larger than
# #' \code{tol_num/minrt}. Default is 50.
#' @param transformtype A character parameter indicating what kind of
#' tranformation('arsinh' or 'none')would be applied to the original matrix.
#' Default is 'none'.
#' @param KNN A interger indicates the number of neighbors in building the
#' KNN graph. Default is 20.
#' @param basecluster A character parameter indicating what kind of simple
#' clustering methods would be used as a preliminary processinng method.
#' Available choices include \code{louvain}, \code{single}, \code{complete},
#' \code{average} and \code{kmeans}. Default is 'louvain'.
#' @param dobasecluster A Bloolean parameter indicates wether to do clustering
#' on the first level. Default is FALSE.
#' @param baseclusternum A integer indicates the number of clusters on the
#' first level. This parameter is useful only when the parameter
#''dobasecluster' is TRUE, and the clustering method is one of  'single',
#''complete', 'average' and 'kmeans'. For 'louvain', we return the division
#' of communities which there are just less communities than this number on
#' the level in the output of function'cluster_louvain' in igraph package.
#' Default is a tenth of the number of samples.
#' @param emphasizedistance A integer indicates the power of element in
#' the C matrix, bigger the emphasizdistance, more likely to detect
#' elongated clusters. Default is 1. More details please refer to the article.
#' @param weightDist A numeric parameter indicates the power of distance.
#' Default is 2.
#' @param weightDens A numeric parameter indicates the power of density.
#' Default is 0.5.
#' @param cutpoint A numeric value in [0,1] indicates the threshold of
#' the quantile of pseudo energy, default is 0.05.
#' @param showprocess A Boolean parameter indicating whether to print
#' intermediate information.Default is FALSE.
#' @param bn A numeric parameter indicates the power of the Minkowski distance.
#' Default is 2.
# #' @param stop_rate A numeric parameter indicating a stopping criterion that
# #' if the number of samples belonging to qualified clusters exceeds
# #' \code{stop_rate*num_data}on this level, then the algorithm will stop and
# #' not bulid the higher hierachical structure. Default is 1.
#' @details The data given by \code{origin_matrix} will be clustered by using
#' the Markov Hierarchical Clustering Algorithm, which generates a hierarchical
#' structure based on the metastability of exponentially perturbed Markov chain.
#' More details of the algorithm please refer to the article.
#' @return This function returns a list including the input parameters, some
#' intermedia results and the hierarchical structure. The list includes the
#' following components:
#' @author Zhenyi Wang
#' @export

MarkovHC = function(origin_matrix,
                    #minrt=50,
                    transformtype="none",
                    KNN=20,
                    basecluster="louvain",
                    dobasecluster=FALSE,
                    baseclusternum=NULL,
                    emphasizedistance=1,
                    weightDist=2,
                    weightDens=0.5,
                    cutpoint=0.01,
                    showprocess=FALSE,
                    bn=2,
                    minBasinSize=0.2,
                    #smallBasinSize=3,
                    noiseBasinSize=10
                    ){
  ##step01.check the input parameters------------------------------------------
  #origin_matrix
  if((!is.matrix(origin_matrix))&(!is.data.frame(origin_matrix))){
    print("The type of 'origin_matrix' should be matirx or dataframe!")
    return(NULL)
  }
  origin_matrix <- t(origin_matrix)
  #minrt
  #if(!is.numeric(minrt)){
  #  print("The type of 'minrt' should be numeric!")
  #  return(NULL)
  #}
  #if(minrt<=2){
  #  print("The value of 'minrt' is too small!")
  #  return(NULL)
  #}

  #transformtype
  if(!is.character(transformtype)){
    print("The type of 'transformtype' should be character!")
    return(NULL)
  }
  if(transformtype=="arsinh"){
    transformed_matrix<-arsinh(origin_matrix)
  }else if(transformtype=="none"){
    transformed_matrix<-origin_matrix
  }else{
    print("This kind of 'transformtype' is unavailable right now!")
    return(NULL)
  }

  #KNN
  if(!is.numeric(KNN)){
    print("The type of 'KNN' should be numeric!")
    return(NULL)
  }

  #basecluster
  if(!is.character(basecluster)){
    print("The type of 'basecluster' should be character!")
    return(NULL)
  }
  if((basecluster=="louvain")|(basecluster=="single")|
     (basecluster=="complete")|(basecluster=="average")|
     (basecluster=="kmeans")){
    basecluster<-basecluster
  }else{
    print("This kind of 'basecluster' is unavailable right now!")
    return(NULL)
  }

  #emphasizedistance
  if(!is.numeric(emphasizedistance)){
    print("The type of 'emphasizedistance' should be numeric!")
    return(NULL)
  }
  if(emphasizedistance<=0){
    print("The parameter 'emphasizedistance' should be a positive number!")
    return(NULL)
  }
  #weightDist
  if(!is.numeric(weightDist)){
    print("The type of 'weightDist' should be numeric!")
    return(NULL)
  }
  #weightDens
  if(!is.numeric(weightDens)){
    print("The type of 'weightDens' should be numeric!")
    return(NULL)
  }
  #cutpoint
  if(!is.numeric(cutpoint)){
    print("The type of 'cutpoint' should be numeric!")
    return(NULL)
  }
  #bn
  if(!is.numeric(bn)){
    print("The type of 'bn' should be numeric!")
    return(NULL)
  }
  #stop_rate
  #if(!is.numeric(stop_rate)){
  #  print("The type of 'stop_rate' should be numeric!")
  #  return(NULL)
  #}
  set.seed(1)
  #Do parallel
  ncore<-detectCores()
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)

  ##step02.calculate SNN and build a strongly connected KNN graph--------------
  dm <- dist(transformed_matrix,method = "minkowski",p=bn)
  dm_matrix <- as.matrix(dm)
  #step02.1 use sNN function to calculate shared nearest neighbors
  #Only calculate sNN of points are in each others snn list.
  sNN_res <- sNN(x=dm,k=KNN,sort = FALSE)

  #step02.2 build a strongly connected KNN graph

  KNN_graph <- matrix(0,nrow(transformed_matrix),nrow(transformed_matrix))
  for (i in 1:nrow(sNN_res$shared)) {
     KNN_graph[i,sNN_res$id[i,]] <- sNN_res$shared[i,]
  }

  #If some people install the develop vesion of dbscan, we should use the
  #block below, refer to https://github.com/mhahsler/dbscan/issues/37 that
  #I reported to them.
  # for(i in 1:nrow(sNN_res$shared)){
  #   KNN_graph[setdiff(1:nrow(KNN_graph),sNN_res$id[i,]),i] <- 0
  # }

  # deprecated, because it's too slow
  # for (i in 1:nrow(dm_matrix)) {
  #   for (j in 1:ncol(dm_matrix)) {
  #     KNN_i <- order(dm_matrix[i,], decreasing = FALSE)[1:KNN]
  #     KNN_j <- order(dm_matrix[j,], decreasing = FALSE)[1:KNN]
  #     KNN_graph[i,j] <- intersect(KNN_i,KNN_j)%>%length()
  #   }
  # }

  # deprecated, because it's not fast enough
  # SNNCalculate <- function(x,y){
  #   intersect(x,y)%>%length()
  # }
  # test_snn <- matrix(0,nrow(transformed_matrix),nrow(transformed_matrix))
  # for (i in 1:nrow(test)) {
  #   test_snn[i,] <- apply(test, 1, SNNCalculate, test[i,])
  # }

  #deprecated, accelerate by Rcpp, but the connections are too strong
  # orderKNN <- function(x,y){
  #   order(x,decreasing = FALSE)[1:y]
  # }
  # KNN_index_temp <- apply(dm_matrix, 1, orderKNN, KNN)
  # KNN_index_temp <- t(KNN_index_temp)
  # print('Calculate shared nearest neighbors of points.')
  # cppFunction('
  # IntegerMatrix SNNCalculate(IntegerMatrix KNN_matrix){
  #   int nRow1 = KNN_matrix.rows();
  #   int nRow2 = KNN_matrix.rows();
  #   IntegerMatrix SNN_matrix(nRow1,nRow2);
  #   Function intersect("intersect");
  #   for(int i=0;i < nRow1;i++){
  #      for(int j=i;j < nRow2;j++){
  #       IntegerVector vector1 = KNN_matrix( i, _);
  #       IntegerVector vector2 = KNN_matrix( j, _);
  #       IntegerVector intersectNN = intersect(_["x"]=vector1,_["y"]=vector2);
  #       int SNN_count = intersectNN.size();
  #       SNN_matrix(i,j) = SNN_count;
  #      }
  #   }
  #   return(SNN_matrix);
  #  }')
  # KNN_graph <- SNNCalculate(KNN_index_temp)
  # for (i in 1:nrow(KNN_graph)) {
  #   KNN_graph[i,setdiff(1:ncol(KNN_graph),KNN_index_temp[i,])] <- 0
  # }
  # KNN_graph[lower.tri(KNN_graph)] <- t(KNN_graph)[lower.tri(KNN_graph)]

  #convert asymmetric matrix of KNN to symmetric matrix of KNN
  #But here we use SNN, this step could be deprecated

  KNN_graph_T <- t(KNN_graph)
  KNN_graph_index <- KNN_graph
  KNN_graph_index[which(KNN_graph_index>0)] <- Inf
  KNN_graph_T_index <- KNN_graph_T
  KNN_graph_T_index[which(KNN_graph_T_index>0)] <- Inf
  symmetric_KNN_graph <- KNN_graph+KNN_graph_T

  #the elements on double positive index should be divided by 2, do not need to do
  #that on one or zero positive index.
  symmetric_KNN_graph[is.infinite(KNN_graph_index)&is.infinite(KNN_graph_T_index)] <- symmetric_KNN_graph[is.infinite(KNN_graph_index)&is.infinite(KNN_graph_T_index)]/2

  rm(list = c('KNN_graph','KNN_graph_T', 'KNN_graph_index', 'KNN_graph_T_index'))
  gc(verbose=FALSE)

  ##below method is deprecated because of slow calculateing speed
  # KNN_graph <- as.matrix(dm)
  # bulid_KNN_row <- function(x,n=KNN){
  #   cut <- x[order(x,decreasing = FALSE)][n+1]
  #   x[which(x)>cut] <- Inf
  #   return(x)
  # }
  # KNN_graph <- apply(X=KNN_graph, MARGIN = 1, bulid_KNN_row)
  # #convert asymmetric matrix to symmetric matrix
  # KNN_graph_T <- t(KNN_graph)
  # KNN_graph_index <- KNN_graph
  # KNN_graph_index[is.finite(KNN_graph_index)] <- 1
  # KNN_graph_index[is.infinite(KNN_graph_index)] <- 0
  # KNN_graph_T_index <- KNN_graph_T
  # KNN_graph_T_index[is.finite(KNN_graph_T_index)] <- 1
  # KNN_graph_T_index[is.infinite(KNN_graph_T_index)] <- 0
  # KNN_graph <- KNN_graph*KNN_graph_index
  # KNN_graph[is.na(KNN_graph)] <- 0
  # KNN_graph_T <- KNN_graph_T*KNN_graph_T_index
  # KNN_graph_T[is.na(KNN_graph_T)] <- 0
  # symmetric_KNN_graph <- KNN_graph+KNN_graph_T
  # symmetric_KNN_graph[(KNN_graph_index==1)&(KNN_graph_T_index==1)] <- symmetric_KNN_graph[(KNN_graph_index==1)&(KNN_graph_T_index==1)]/2
  # symmetric_KNN_graph[symmetric_KNN_graph==0] <- Inf
  # diag(symmetric_KNN_graph) <- 0

  ##other options are calculating the degree or PageRank of each node in the graph, we take the degree as the density of the node here.
  #use 'degree' in igraph to calculate the degree of vertexes in the network.

  #symmetric_KNN_graph_sparse <- as(symmetric_KNN_graph, "dgCMatrix") %>% Matrix::summary() %>% as.data.frame()
  #symmetric_KNN_graph_object <- make_graph(t(symmetric_KNN_graph_sparse[,1:2]), directed = FALSE)
  #graph_attr(symmetric_KNN_graph_object,'weight') <- symmetric_KNN_graph_sparse[,3]
  #use graph_from_adjacency_matrix for convenience
  symmetric_KNN_graph_object <- graph_from_adjacency_matrix(adjmatrix = symmetric_KNN_graph,
                                                            mode = 'undirected',
                                                            weighted = TRUE,
                                                            diag = TRUE)

  #centrality_scores <- eigen_centrality(symmetric_KNN_graph_object, weights = symmetric_KNN_graph_sparse[,3])$vector
  centrality_scores <- degree(symmetric_KNN_graph_object, v = V(symmetric_KNN_graph_object),
                              mode = "total",
                              loops = TRUE, normalized = FALSE)+1

  ##step03.do preclustering----------------------------------------------------
  #hierarchical clustering or k-means clutering or finding Maximum clique.
  if(dobasecluster==TRUE){
    #do clustering on the first level
    #Use one type of hierarchical clustering as the basic clustering tool
    if((basecluster=="single")|(basecluster=="complete")|(basecluster=="average")){
      hresult <- hclust(dm,method = basecluster)
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
      }
        #user sets the baseclusternum parameter
        hresult_cut <- cutree(hresult,k=baseclusternum)
    }else if(basecluster=="kmeans"){
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
      }
        #user sets the baseclusternum parameter
        kmeansresult <- kmeans(transformed_matrix, centers=baseclusternum,iter.max = 500)
        hresult_cut <- kmeansresult$cluster
    }else if(basecluster=='louvain'){
        # #find all maximum cliques on the graph
        # maxcliques <- max_cliques(symmetric_KNN_graph_object)
        # #regard a clique as a cluster
        # hresult_cut <- integer(length = nrow(transformed_matrix))
        # for (index_maxcliques in 1:length(maxcliques)) {
        #   hresult_cut[as.integer(maxcliques[[index_maxcliques]])] <- index_maxcliques
        # }
        if (is.null(baseclusternum)){
          baseclusternum <- ceiling(nrow(transformed_matrix)/10)
        }
        cluster_louvain_object <- cluster_louvain(graph = symmetric_KNN_graph_object,
                                                  weights = E(symmetric_KNN_graph_object)$weight )
        for (i in 1:nrow(cluster_louvain_object$memberships)) {
          if(length(unique(cluster_louvain_object$memberships[i,]))<=baseclusternum){
            hresult_cut <- cluster_louvain_object$memberships[i,]
            break
            }
        }
     }
    }else{
     #do not do clustering on the first level
     hresult_cut <- 1:nrow(transformed_matrix)
    }
    #merge a cluster as a single point
    #or
    #merge a clique as a single point
    #every similarity to other clusters is the maximum similarity between the cluster with other clusters
    unique_clusters <- unique(hresult_cut)
    #merge rows
    symmetric_KNN_graph_merged <- matrix(0,length(unique_clusters),nrow(transformed_matrix))
    for(clusterindex in 1:length(unique_clusters)){
      temp_index <- which(hresult_cut==clusterindex)
      if(length(temp_index)==1){
        symmetric_KNN_graph_merged[clusterindex,] <- symmetric_KNN_graph[temp_index, ]
      }else{
        temp_cluster <- symmetric_KNN_graph[temp_index, ]
        symmetric_KNN_graph_merged[clusterindex,] <- apply(temp_cluster, 2, max)
      }
    }
    #merge columns
    #each elements in symmetric_KNN_graph_cluster is the similarity bwt clusters and clusters
    symmetric_KNN_graph_cluster <- matrix(0,nrow(symmetric_KNN_graph_merged),nrow(symmetric_KNN_graph_merged))
    for(clusterindex in 1:length(unique_clusters)){
      for(clusterindex2 in 1:length(unique_clusters)){
        temp_cluster <- symmetric_KNN_graph_merged[clusterindex, which(hresult_cut==clusterindex2)]
        symmetric_KNN_graph_cluster[clusterindex, clusterindex2] <- max(temp_cluster)
      }
    }
    #diag(symmetric_KNN_graph_cluster) <- min(symmetric_KNN_graph_cluster)
    diag(symmetric_KNN_graph_cluster) <- KNN
    #calculate the centrality_scores of clusters
    centrality_scores_cluster <- integer(length = length(unique_clusters))
    for (score_index in 1:length(unique_clusters)) {
      centrality_scores_tpm <- centrality_scores[which(hresult_cut==score_index)]
      centrality_scores_cluster[score_index] <- mean(centrality_scores_tpm)
    }

  ## Main part of MarkovHC algorithm
  ##step04. Calculate the transition probability matrix and the pseudo energy matrix
  #step04.1 Calculate the transition probability matrix
  transitionMatrix<-transition_probability(matrix=symmetric_KNN_graph_cluster,
                                           densevector=centrality_scores_cluster,
                                           weightDist=weightDist,
                                           weightDens=weightDens)

  #step04.2 Calculate the pseudo energy matrix
  C_matrix <- Calculate_C_Matrix(matrix=symmetric_KNN_graph_cluster,
                                 densevector=centrality_scores_cluster,
                                 emphasizedistance=emphasizedistance,
                                 weightDist=weightDist,
                                 weightDens=weightDens)
  C_matrix <- C_matrix + 0.1
  C_matrix[which(C_matrix==Inf)] <- 0
  # C_matrix_graph_sparse <- as(as.matrix(C_matrix), "dgCMatrix")%>%Matrix::summary()%>%as.data.frame()
  # C_matrix_graph_object <- make_graph(t(C_matrix_graph_sparse[,1:2]), directed = TRUE)
  # graph_attr(C_matrix_graph_object,'weight') <- C_matrix_graph_sparse[,3]
  C_matrix_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(C_matrix),
                                                       mode = 'directed',
                                                       weighted = TRUE,
                                                       diag = TRUE)
  #step04.3 Calculate the shortest distance between each vertex pair in the graph
  print('Calculate the shortest distance between each vertex pair in the graph.')
  C_matrix_graph_shortest_distance <- matrix(data = 0, nrow = nrow(C_matrix), ncol = nrow(C_matrix))
  for (i in 1:nrow(C_matrix)) {
    for (j in 1:nrow(C_matrix)) {
      if(i==j){next}
      C_matrix_graph_shortest_distance[i,j] <- distances(graph = C_matrix_graph_object,
                                                         v = i,
                                                         to = j,
                                                         mode = 'out',
                                                         weights = E(C_matrix_graph_object)$weight,
                                                         algorithm = "dijkstra")%>%as.vector()
    }
  }

  ##step05. Build the hierarchical structure-----------------------------------
  P_updated <- transitionMatrix
  MarkovHC_result <- list()
  #Store the result of base clustering
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
                       basinNum=basinNum)
  MarkovHC_result <- append(MarkovHC_result, list(level_result))
  levels_indice <- 1
  print(paste('Build the level ',as.character(levels_indice),'...', sep = ''))
  while (TRUE) {
    levels_indice <- levels_indice + 1
    print(paste('Build the level ',as.character(levels_indice),'...', sep = ''))
    ##step05.1 find basins and attractors
    RS_vector <- judge_RS(P=P_updated)

    ##step05.2 constructe the list to store the result of this level
    attractors <- list()
    basins <- list()
    graphvertex_attractors <- list()
    graphvertex_basins <- list()
    attractorPoints <- list()
    basinPoints <- list()

    ##step05.3 partition the state space
    processed_attractors <- integer(length = length(RS_vector))
    #processed_attractors[which(RS_vector==0)] <- 1
    processed_attractors[which(RS_vector<=1e-5)] <- 1
    attractor_indice <- 1
    while(TRUE){
     if(all(processed_attractors==1)){break}
     print(paste('Find attractors in the basin ',as.character(attractor_indice),'.', sep = ''))
     attractor_temp <- which(processed_attractors==0)[1]
     #processed_attractors[attractor_temp] <- 1
     # P_updated_graph_sparse <- as(as.matrix(P_updated), "dgCMatrix")%>%Matrix::summary()%>%as.data.frame()
     # P_updated_graph_object <- make_graph(t(P_updated_graph_sparse[,1:2]), directed = TRUE)
     # graph_attr(P_updated_graph_object,'weight') <- P_updated_graph_sparse[,3]
     P_updated_graph_object <- graph_from_adjacency_matrix(adjmatrix = as.matrix(P_updated),
                                                           mode = 'directed',
                                                           weighted = TRUE,
                                                           diag = TRUE)
     #attractor_temp_access <- all_simple_paths(graph = P_updated_graph_object, from = attractor_temp,
     #                                         mode = 'out')%>%unlist()%>%unique()
     #subcomponent includes the vertex itself in this graph
     attractor_temp_access <- subcomponent(graph = P_updated_graph_object,
                                           v = attractor_temp,
                                           mode = "out")%>%as.vector()
     processed_attractors[attractor_temp_access] <- 1
     #attractors <- c(attractors, list(unique(c(attractor_temp_access, attractor_temp))))
     attractors <- c(attractors, list(attractor_temp_access))
     #basins_temp_merged <- c()
     #  for (i in unique(c(attractor_temp_access, attractor_temp))) {
     #    basins_temp <- all_simple_paths(graph = P_updated_graph_object, from = i,
     #                                    mode = 'in')%>%unlist()%>%unique()
     #   basins_temp_merged <- c(basins_temp_merged, basins_temp)
     #  }
     #basins_temp_merged <- c(basins_temp_merged, unique(c(attractor_temp_access, attractor_temp)))%>%unique()
     basins_temp_merged <- subcomponent(graph = P_updated_graph_object,
                                        v = attractor_temp,
                                        mode = "in")%>%as.vector()
     basins <- c(basins, list(basins_temp_merged))
     attractor_indice <- attractor_indice+1
    }

    ##step05.4 assign the points and graph vertexes to basins and attractors
    basinNum <- length(basins)
    basin_indice <- 1
    for (i in 1:basinNum){
        print(paste('Partition the basin ',as.character(basin_indice),'.', sep = ''))
        #assign attractor points
        attractorPoints_temp <- level_result$attractorPoints[attractors[[i]]]%>%unlist()%>%unique()
        attractorPoints <- c(attractorPoints, list(attractorPoints_temp))
        #assign basin points
        basinPoints_temp <- level_result$basinPoints[basins[[i]]]%>%unlist()%>%unique()
        basinPoints <- c(basinPoints, list(basinPoints_temp))
        #assign graphvertex_attractors
        graphvertex_attractors_temp <- level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique()
        graphvertex_attractors <- c(graphvertex_attractors, list(graphvertex_attractors_temp))
        #assign graphvertex_basins
        graphvertex_basins_temp <- level_result$graphvertex_basins[basins[[i]]]%>%unlist()%>%unique()
        graphvertex_basins <- c(graphvertex_basins, list(graphvertex_basins_temp))
        basin_indice <- basin_indice+1
    }

    ##step05.5 update the pseudo energy matrix
    print('Update the pseudo energy matrix.')
    C_matrix_updated <- matrix(data = 0, nrow = basinNum, ncol = basinNum)
    #for (i in 1:basinNum) {
    #  for (j in 1:basinNum) {
    #    if(i==j){next}
    #    C_matrix_updated[i,j] <- distances(graph = C_matrix_graph_object,
    #                                       v = level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique(),
    #                                       to = level_result$graphvertex_basins[basins[[j]]]%>%unlist()%>%unique(),
    #                                       mode = 'out',
    #                                       weights = E(C_matrix_graph_object)$weight,
    #                                       algorithm = "dijkstra")%>%min()
    #  }
    #}
    for (i in 1:basinNum) {
     for (j in 1:basinNum) {
       if(i==j){next}
       source_v <- level_result$graphvertex_attractors[attractors[[i]]]%>%unlist()%>%unique()
       target_v <- level_result$graphvertex_basins[basins[[j]]]%>%unlist()%>%unique()
       C_matrix_updated[i,j] <- C_matrix_graph_shortest_distance[source_v, target_v]%>%min()
     }
    }

    #If the basin is small, we regard it as noise and merge it to the closest basin
    C_matrix_updated_mergedsmallbasin <- C_matrix_updated

    noise_basins <- c()
    max_basin_size <- c()
    #find noise basins and record the size of the biggest basin
    for (max_basin_indice in 1:length(basinPoints)) {
      max_basin_size <- max(max_basin_size, length(basinPoints[[max_basin_indice]]))
      if(length(basinPoints[[max_basin_indice]])<=noiseBasinSize){
        noise_basins <- c(noise_basins, max_basin_indice)
      }
    }

    #qualified basins
    qualified_basins <- setdiff(1:nrow(C_matrix_updated),noise_basins)
    #calculate distances bwt noise basins with qualified basins
    noise_basins_to_qualified_basins <- matrix(Inf,nrow(C_matrix_updated),nrow(C_matrix_updated))
    for (noise_basins_i in noise_basins) {
      for (qualified_basins_i in qualified_basins) {
        noise_basins_to_qualified_basins_temp <- dm_matrix[basinPoints[[noise_basins_i]],basinPoints[[qualified_basins_i]]]
        noise_basins_to_qualified_basins[noise_basins_i,qualified_basins_i] <- min(noise_basins_to_qualified_basins_temp)
      }
    }

    #modify energy matrix
    if((max_basin_size/nrow(transformed_matrix))>=minBasinSize){
      print('Merge noise basins to qualified basins.')
      row_min <- apply(C_matrix_updated_mergedsmallbasin, 1, min)
      #let them cannot be recurrent
      for (noise_basins_i in noise_basins) {
        C_matrix_updated_mergedsmallbasin[noise_basins_i,noise_basins_i] <- Inf
      }
      #merge noise basins
      for (noise_basins_i in noise_basins) {
        #this noise_basin is an outlier which cannot reach any other basins
        #merge it to quilified basins
        if(sum(is.infinite(C_matrix_updated_mergedsmallbasin[noise_basins_i,]))==(nrow(C_matrix_updated))){
          cloest_basin_index <- which(noise_basins_to_qualified_basins[noise_basins_i,]==min(noise_basins_to_qualified_basins[noise_basins_i,]))
          C_matrix_updated_mergedsmallbasin[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
        }else{
          #this noise_basin is not an outlier
          #let it merge its cloest basins
          cloest_basin_index <- which(C_matrix_updated_mergedsmallbasin[noise_basins_i,]==min(C_matrix_updated_mergedsmallbasin[noise_basins_i,]))
          C_matrix_updated_mergedsmallbasin[noise_basins_i,cloest_basin_index] <- row_min[noise_basins_i]
        }
      }
    }

    ##step05.6 update the transition probability matrix
    C_matrix_updated_indice <- C_matrix_updated_mergedsmallbasin
    diag(C_matrix_updated_indice) <- Inf
    if((basinNum>1)&(all(is.infinite(C_matrix_updated_indice))==FALSE)){
      print('Update the transition probability matrix.')
      update_P_result <- update_P(C_matrix_updated=C_matrix_updated_mergedsmallbasin, C_cut=cutpoint)
      P_updated <- update_P_result[[1]]
      energyCutpoint <- update_P_result[[2]]
    }

    ##step05.7 constructe the list to store the result of MarkovHC algorithm
    level_result <- list(basins=basins,
                         attractors=attractors,
                         graphvertex_attractors=graphvertex_attractors,
                         graphvertex_basins=graphvertex_basins,
                         basinPoints=basinPoints,
                         attractorPoints=attractorPoints,
                         basinNum=basinNum,
                         C_matrix_updatedmatrix=C_matrix_updated,
                         C_matrix_updated_mergedsmallbasin=C_matrix_updated_mergedsmallbasin,
                         P_updated=P_updated,
                         energyCutpoint=energyCutpoint
                         )
    MarkovHC_result <- append(MarkovHC_result, list(level_result))

    if((basinNum==1)|(all(is.infinite(C_matrix_updated_indice)))){
      ##step06. Output the results---------------------------------------------
      #The input parameters
      inputParameters <- list(
        #minrt=minrt,
        transformtype=transformtype,
        KNN=KNN,
        basecluster=basecluster,
        dobasecluster=dobasecluster,
        baseclusternum=baseclusternum,
        emphasizedistance=emphasizedistance,
        weightDist=weightDist,
        weightDens=weightDens,
        cutpoint=cutpoint,
        showprocess=showprocess,
        bn=bn
        #stop_rate=stop_rate
      )
      #The results among the process
      midResults <- list(
        symmetric_KNN_graph = symmetric_KNN_graph,
        symmetric_KNN_graph_object = symmetric_KNN_graph_object,
        sNN_res = sNN_res,
        centrality_scores = centrality_scores,
        symmetric_KNN_graph_cluster = symmetric_KNN_graph_cluster,
        transitionMatrix = transitionMatrix,
        dm_matrix = dm_matrix,
        C_matrix = C_matrix,
        C_matrix_graph_object=C_matrix_graph_object,
        C_matrix_graph_shortest_distance = C_matrix_graph_shortest_distance,
        centrality_scores_cluster = centrality_scores_cluster
      )
      #The MarkovHC object
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
