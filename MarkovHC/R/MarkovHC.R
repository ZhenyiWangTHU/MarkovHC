#' Markov Hierarchical Clustering Algorithm
#'
#' Function \code{MarkovHC} The main function of MarkovHC, which applies Markov
#' Hierarchical Clustering Algorithm to the given data (a dataframe or a matrix).
#' @param origin_matrix A matrix or a dataframe with row names and column names.
#' Each row should be a feature and each column should be a sample.
#' @param minrt A parameter helping to judge whether a group of samples is a
#' qualified cluster. If we assume that we have \code{tol_num} samples in
#' total, then a qualified cluster should be of size larger than
#' \code{tol_num/minrt}. Default is 50.
#' @param transformtype A character parameter indicating what kind of
#' tranformation('arsinh' or 'none')would be applied to the original matrix.
#' Default is 'none'.
#' @param KNN A interger indicates the number of neighbors in building the
#' KNN graph. Default is 20.
#' @param basecluster A character parameter indicating what kind of simple
#' clustering methods would be used as a preliminary processinng method.
#' Available choices include \code{clique}, \code{single}, \code{complete},
#' \code{average} and \code{kmeans}. Default is 'clique'.
#' @param dobasecluster A Bloolean parameter indicates wether to do clustering
#' on the first level. Default is FALSE.
#' @param baseclusternum A integer indicates the number of clusters on the
#' first level. This parameter is useful only when the parameter
#''dobasecluster' is TRUE. Default is a tenth of the number of samples.
#' @param emphasizedistance A integer indicates the power of element in
#' the C matrix, bigger the emphasizdistance, more likely to detect
#' elongated clusters. Default is 1. More details please refer to the article.
#' @param weightDist A numeric parameter indicates the power of distance.
#' Default is 2.
#' @param weightDens A numeric parameter indicates the power of density.
#' Default is 0.5.
#' @param showprocess A Boolean parameter indicating whether to print
#' intermediate information.Default is FALSE.
#' @param bn A numeric parameter indicates the power of the Minkowski distance.
#' Default is 2.
#' @param stop_rate A numeric parameter indicating a stopping criterion that
#' if the number of samples belonging to qualified clusters exceeds
#' \code{stop_rate*num_data}on this level, then the algorithm will stop and
#' not bulid the higher hierachical structure. Default is 1.
#' @details The data given by \code{origin_matrix} will be clustered by using
#' the Markov Hierarchical Clustering Algorithm, which generates a hierarchical
#' structure based on the metastability of exponentially perturbed Markov chain.
#' More details of the algorithm please refer to the article.
#' @return This function returns a list including the input parameters, some
#' intermedia results and the hierarchical structure. The list includes the
#' following components:
#' @author Zhenyi Wang
#' @export

MarkovHC<-function(origin_matrix,
              minrt=50,
              transformtype="none",
              KNN=20,
              basecluster="clique",
              dobasecluster=FALSE,
              baseclusternum=NULL,
              emphasizedistance=1,
              weightDist=2,
              weightDens=0.5,
              showprocess=FALSE,
              bn=2,
              stop_rate=1){
  ##step01.check the input parameters------------------------------------------
  #origin_matrix
  if((!is.matrix(origin_matrix))&(!is.data.frame(origin_matrix))){
    print("The type of 'origin_matrix' should be matirx or dataframe!")
    return(NULL)
  }
  origin_matrix <- t(origin_matrix)
  #minrt
  if(!is.numeric(minrt)){
    print("The type of 'minrt' should be numeric!")
    return(NULL)
  }
  if(minrt<=2){
    print("The value of 'minrt' is too small!")
    return(NULL)
  }

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
  if((basecluster=="clique")|(basecluster=="single")|
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
  #bn
  if(!is.numeric(bn)){
    print("The type of 'bn' should be numeric!")
    return(NULL)
  }
  #stop_rate
  if(!is.numeric(stop_rate)){
    print("The type of 'stop_rate' should be numeric!")
    return(NULL)
  }

  #Do parallel
  ncore<-detectCores()
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)

  ##step02.calculate SNN and build a strongly connected KNN graph--------------
  dm <- dist(transformed_matrix,method = "minkowski",p=bn)
  #step02.1 use sNN function to calculate shared nearest neighbors
  sNN_res <- sNN(x=dm,k=KNN,sort = FALSE)
  #step02.2 build a strongly connected KNN graph
  KNN_graph <- matrix(0,nrow(dm),ncol(dm))
  for (i in 1:nrow(sNN_res$shared)) {
    KNN_graph[i,sNN_res$id[i,]] <- sNN_res$shared[i,]
  }
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
  #use 'eigen_centrality' in igraph to find Eigenvector Centrality Scores of Network Positions.
  symmetric_KNN_graph_sparse <- as(as.matrix(symmetric_KNN_graph), "dgCMatrix")%>%summary()%>%as.data.frame()
  symmetric_KNN_graph_object <- make_graph(t(symmetric_KNN_graph_sparse[,1:2]), directed = FALSE)
  graph_attr(symmetric_KNN_graph_object,'weight') <- symmetric_KNN_graph_sparse[,3]
  centrality_scores <- eigen_centrality(symmetric_KNN_graph_object, weights = symmetric_KNN_graph_sparse[,3])$vector

  ##step03.do preclustering----------------------------------------------------
  #hierarchical clustering or k-means clutering or finding Maximum clique.
  if(dobasecluster==TRUE){
    #do clustering on the first level
    #Use one type of hierarchical clustering as the basic clustering tool
    if((basecluster=="single")|(basecluster=="complete")|(basecluster=="average")){
      hresult<-hclust(dm,method = basecluster)
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
        hresult_cut <- cutree(hresult,k=baseclusternum)
      }else{
        #user sets the baseclusternum parameter
        hresult_cut <- cutree(hresult,k=baseclusternum)
      }
    }else if(basecluster=="kmeans"){
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
        kmeansresult <- kmeans(transformed_matrix, centers=baseclusternum,iter.max = 500)
        hresult_cut <- kmeansresult$cluster
      }else{
        #user sets the baseclusternum parameter
        kmeansresult <- kmeans(transformed_matrix, centers=baseclusternum,iter.max = 500)
        hresult_cut <- kmeansresult$cluster
      }
    }else if(basecluster=='clique'){
        #find all maximum cliques on the graph
        maxcliques <- max_cliques(symmetric_KNN_graph_object)
        #regard a clique as a cluster
        hresult_cut <- integer(length = nrow(transformed_matrix))
        for (index_maxcliques in 1:length(maxcliques)) {
          hresult_cut[as.integer(maxcliques[[index_maxcliques]])] <- index_maxcliques
        }
     }
    }else{
     #do not do clustering on the first level
     hresult_cut <- 1:nrow(transformed_matrix)
    }
    #merge a cluster as a single point
    #or
    #merge a clique as a single point
    #every distance to out clusters is the minimum distance from the cluster to out clusters
    unique_clusters <- unique(hresult_cut)
    #replace the 0 in graph with Inf, convert the graph matrix to the similarity matrix
    symmetric_KNN_graph_similarity <- symmetric_KNN_graph
    symmetric_KNN_graph_similarity[which(symmetric_KNN_graph_similarity)==0] <- Inf
    #merge rows
    symmetric_KNN_graph_similarity_merged <- matrix(0,length(unique_clusters),nrow(transformed_matrix))
    for(clusterindex in 1:length(unique_clusters)){
      temp_cluster <- symmetric_KNN_graph_similarity[which(hresult_cut==clusterindex), ]
      symmetric_KNN_graph_similarity_merged[clusterindex,] <- apply(temp_cluster, 2, min)
    }
    #merge columns
    #each elements in symmetric_KNN_graph_similarity_cluster is the distance bwt clusters and clusters
    symmetric_KNN_graph_similarity_cluster <- matrix(0,nrow(symmetric_KNN_graph_similarity_merged),nrow(symmetric_KNN_graph_similarity_merged))
    for(clusterindex in 1:length(unique_clusters)){
      for(clusterindex2 in 1:length(unique_clusters)){
        temp_cluster <- symmetric_KNN_graph_similarity_merged[clusterindex, which(hresult_cut==clusterindex2)]
        symmetric_KNN_graph_similarity_cluster[clusterindex, clusterindex2] <- min(temp_cluster)
      }
    }

    #calculate the centrality_scores of clusters
    centrality_scores_cluster <- integer(length = length(unique_clusters))
    for (score_index in 1:length(unique_clusters)) {
      centrality_scores_tpm <- centrality_scores[which(hresult_cut==score_index)]
      centrality_scores_cluster[score_index] <- mean(centrality_scores_tpm)
    }

  ## Main MarkovHC algorithm
  ##step04. Calculate the transition probability matrix and the pseudo energy matrix
  #step04.1 Calculate the transition probability matrix
  transitionMatrix<-transition_probability(matrix=symmetric_KNN_graph_similarity,
                                           densevector=centrality_scores_cluster,
                                           weightDist=weightDist,
                                           weightDens=weightDens)

  #step04.2 Calculate the pseudo energy matrix
  C_matrix<-Calculate_C_Matrix(matrix=symmetric_KNN_graph_similarity,
                               densevector=centrality_scores_cluster,
                               emphasizedistance=emphasizedistance,
                               weightDist=weightDist,
                               weightDens=weightDens)

  ##step05. Build the hierarchical structure-----------------------------------
  P_updated <- transitionMatrix
  while (TRUE) {
    ##step05.1 find basins and attractors
    RS_vector <- judge_RS(P=transitionMatrix)
    processed_attractors <- integer(length = length(RS_vector))
    while(TRUE){


    }


    ##step05.2 update the pseudo energy matrix



    ##step05.3 update the transition probability matrix
    #find the shortest path https://igraph.org/r/doc/distances.html


    ##step05.4 constructe the list to store the result of this level
    basinsPoints <- list()
    attractorsPoints <- list()
    energyMatrix <- matrix()
    transMatrix <- matrix()
    basinNum <- length(basinsPoints)
    level_result <- list(basins=basinsPoints,
                         attractors=attractorsPoints,
                         energyMatrix=energyMatrix,
                         transMatrix=transMatrix,
                         basinNum=basinNum)
    ##step05.5 constructe the list to store the result of MarkovHC algorithm
    MarkovHC_result <- c(MarkovHC_result, level_result)
  }




}
