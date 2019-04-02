#' Markov Hierarchical Clustering Algorithm
#'
#' Function \code{MarkovHC} The main function of MarkovHC, which applies Markov
#' Hierarchical Clustering Algorithm to the given data (a dataframe or a matrix).
#' @param origin_matrix A matrix or a dataframe. Each row should be a feature
#' and each column should be a sample.
#' @param minrt A parameter helping to judge whether a group of samples is a
#' qualified cluster. If we assume that we have \code{tol_num} samples in
#' total, then a qualified cluster should be of size larger than
#' \code{tol_num/minrt}. Default is 50.
#' @param transformtype A character parameter indicating what kind of
#' tranformation('arsinh' or 'none')would be applied to the original matrix.
#' Default is 'none'.
#' @param basecluster A character parameter indicating what kind of simple
#' clustering methods would be used as a preliminary processinng method.
#' Available choices include \code{single}, \code{complete}, \code{average}
#' and \code{kmeans}. Default is 'average'.
#' @param emphasizedistance A interger indicates the power of element in
#' the C matrix, bigger the emphasizdistance, more likely to detect
#' elongated clusters. Default is 1. More details please refer to the article.
#' @param weightDist A numeric parameter indicates the power of distance.
#' Default is 2.
#' @param weightDens A numeric parameter indicates the power of density.
#' Default is 0.5.
#' @param showprocess A Boolean parameter indicating whether to print
#' intermediate information.
#' @param bn A numeric parameter indicates the power of the Minkowski distance.
#' Default is 2.
#' @param stop_rate A numeric parameter indicating a stopping criterion that
#' if the number of samples belonging to qualified clusters exceeds
#' \code{stop_rate*num_data}on this level, then the algorithm will stop and
#' not bulid the higher hierachical structure.
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
              basecluster="average",
              emphasizedistance=1,
              weightDist=2,
              weightDens=0.5,
              showprocess=FALSE,
              bn=2,
              stop_rate=0.9){
  #check the input parameters--------------------------------------------------
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

  #basecluster
  if(!is.character(basecluster)){
    print("The type of 'basecluster' should be character!")
    return(NULL)
  }
  if((basecluster=="single")|(basecluster=="complete")|
     (basecluster=="average")|(basecluster=="kmeans")){
    basecluster<-basecluster
  }else{
    print("This kind of 'basecluster' is unavaliable right now!")
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
##20180402
  #Do parallel-----------------------------------------------------------------
  ncore<-detectCores()
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)
  #memory.limit(3000)
  #open parallel calculation
  if(transformtype=="arsinh"){
    transformed_datatable<-arsinh(origin_matrix)
  }else if(transformtype=="none"){
    transformed_datatable<-origin_matrix
  }else{
    print("No such kind of transformation provided.")
    return(NULL)
  }

  ##Preclustering---------------------------------------------------------------
  #Use one type of hierarchical clustering as the basic clustering tool
  if((basecluster=="single")|(basecluster=="complete")|(basecluster=="average")){
    hresult<-hclust(dist(transformed_datatable,method = "minkowski",p=bn),method = basecluster)
    if (nrow(transformed_datatable)>1000){
      hresult_cut <- cutree(hresult,k=1000)
    }else{
      #Small dataset no need base clustering
      hresult_cut <- cutree(hresult,k=nrow(transformed_datatable))
    }
  }else{
    if (nrow(transformed_datatable)>1000){
      kmeansresult <- kmeans(transformed_datatable, centers=1000,iter.max = 500)
      hresult_cut <- kmeansresult$cluster
    }else{
      #Small dataset no need base clustering
      hresult_cut <- 1:nrow(transformed_datatable)
    }
  }
  #down sample data and retain the shape of each cluster
  mingdan <- get_cluster_index(clusterIndex=hresult_cut)
  statematrix <- get_statematrix(matrix=transformed_datatable, basecluster=hresult_cut)

  rm_result<-get_scatter(mingdan,transformed_datatable,statematrix)
  statematrix<-rm_result[[1]]
  mingdan<-rm_result[[2]]
  record_cut<-rm_result[[3]]

  #Calculate the density of each state-----------------------------------------
  if(nrow(statematrix)<1000){
    densevector<-get_densityvector(matrix = statematrix)
  }else{
    densevector<-get_densityvector2(matrix1=statematrix,matrix2=transformed_datatable)
  }

  ##Build the hierarchical structure-------------------------------------------
  #Find recurrent classes on the first level-----------------------------------
  dm <- as.matrix(dist(statematrix,method = "minkowski",p=bn))
  transitionMatrixOnFirstLevel<-transition_first_level(matrix=dm,
                                                       densevector=densevector,
                                                       weightDist=weightDist,
                                                       weightDens=weightDens)
  transitionMatrixOnFirstLevel<-modify_P(transitionMatrixOnFirstLevel,record_cut,densevector)
  mediandm<-median(dm)

  #Calculate C matrix----------------------------------------------------------
  C<-Calculate_C_Matrix(GD=dm,
                        densevector=densevector,
                        emphasizedistance=emphasizedistance,
                        weightDist=weightDist,
                        weightDens=weightDens)
  C<-modify_C(C,record_cut,densevector)

  #Find attractors and basins on the first level-------------------------------
  resultOnTheFirstLevel<-rdv(P=transitionMatrixOnFirstLevel,
                             level=1,
                             densevector=densevector,
                             showprocess=showprocess)
  allresult<-list()
  allresult[[1]]<-resultOnTheFirstLevel

  #Find the shortest path
  ips<-Inner_Path_Dijk(output=resultOnTheFirstLevel,
                       C=C)
  #The inner path
  ip<-ips[["ip"]]
  newlist<-ips[["newlist"]]
  #Find the target boundary
  TBresult<-find_TB(ip,C,-1)
  #Generate the transition matrix on level1
  C_jump <- c()
  GP_result<-generateTransitionMatirx1(TBresult=TBresult,
                                       output=resultOnTheFirstLevel)
  TBresult<-find_TB(ip,C,GP_result[[2]])
  GP_result<-generateTransitionMatirx1(TBresult=TBresult,
                                       output=resultOnTheFirstLevel)
  P <- GP_result[[1]]
  C_jump <- c(C_jump, GP_result[[2]])
  C_distribution <- list(GP_result[[3]])
  #The number of basins on level1
  current_m<-resultOnTheFirstLevel[[1]]$number
  #minimum ratio to be a cluster
  clusterQualification<-ceiling(nrow(statematrix)/minrt)
  #the Current level
  current_level <- 1
  if(showprocess){
    print(paste("Building the level", as.character(current_level), sep = " "))
  }
  #The biggest number of basins on a level in this structure
  bigmark<-0
  #The level index
  bj<-2
  #This loop aims to find attractors and basins on higher levels and bulid the-
  #whole hierarchical structure------------------------------------------------
  while (TRUE){
    #如果合格的大类只剩下一个，而且层数超过5了，而且巅峰大类数超过2，而且被合进合格大类的点数超过了stop_rate这个比例
    if(length(which(lengths(resultOnTheFirstLevel[[2]])>=clusterQualification))<=1&bj>5&bigmark>2&(length(unique(unlist(resultOnTheFirstLevel[[2]][which(lengths(resultOnTheFirstLevel[[2]])>=clusterQualification)])))>(stop_rate*nrow(statematrix)))){break}
    current_level <- current_level+1
    if(showprocess){
      print(paste("Building the level", as.character(current_level), sep = " "))
    }
    outputt2<-buildHigherLevel(output=resultOnTheFirstLevel,
                               P=P,
                               showprocess=showprocess)

    if(resultOnTheFirstLevel[[1]]$number==outputt2[[1]]$number){
      GP_result<-generate_P(soresult=TBresult,
                            output=resultOnTheFirstLevel,
                            oldoutput=allresult[[bj-2]],
                            totalnumber=nrow(statematrix),
                            strengthen=TRUE,
                            showprocess=showprocess)
      P <- GP_result[[1]]
      outputt2<-buildHigherLevel(resultOnTheFirstLevel,P)
      if(resultOnTheFirstLevel[[1]]$number==outputt2[[1]]$number){
        GP_result<-generate_P(soresult=TBresult,
                              output=resultOnTheFirstLevel,
                              oldoutput=allresult[[bj-2]],
                              totalnumber=nrow(statematrix),
                              strengthen=TRUE,
                              mark=1,
                              showprocess=showprocess)
        P <- GP_result[[1]]
        outputt2<-buildHigherLevel(resultOnTheFirstLevel,P)
      }
    }

    newresult<-find_shortest_out_dijk_plus(C=C,
                                           output=outputt2,
                                           oldoutput=resultOnTheFirstLevel,
                                           oldTBresult=TBresult,
                                           oldlist=newlist,
                                           cp=-1)
    TBresult_temp<-newresult[[1]]
    newlist_temp<-newresult[[2]]
    resultOnTheFirstLevel_temp<-outputt2
    allresult[[bj]]<-resultOnTheFirstLevel_temp
    bj<-bj+1
    if(showprocess){
      message(paste("from level",as.character(bj-1),"to level",as.character(bj),sep = " "))
    }
    GP_result<-generate_P(soresult=TBresult_temp,
                          output=resultOnTheFirstLevel_temp,
                          oldoutput=allresult[[bj-2]],
                          totalnumber=nrow(statematrix),
                          showprocess=showprocess)

    newresult<-find_shortest_out_dijk_plus(C=C,
                                           output=outputt2,
                                           oldoutput=resultOnTheFirstLevel,
                                           oldTBresult=TBresult,
                                           oldlist=newlist,
                                           GP_result[[2]])
    TBresult<-newresult[[1]]
    newlist<-newresult[[2]]
    resultOnTheFirstLevel<-outputt2
    GP_result<-generate_P(soresult=TBresult,
                          output=resultOnTheFirstLevel,
                          oldoutput=allresult[[bj-2]],
                          totalnumber=nrow(statematrix),
                          showprocess=showprocess)

    P <- GP_result[[1]]
    C_jump <- c(C_jump, GP_result[[2]])
    C_distribution <- c(C_distribution, list(GP_result[[3]]))
    current_m<-resultOnTheFirstLevel[[1]]$number
    if(current_m==1){
      break
    }
    #update bigmark
    if(length(which(lengths(resultOnTheFirstLevel[[2]])>=clusterQualification))>bigmark){
      bigmark<-length(which(lengths(resultOnTheFirstLevel[[2]])>=clusterQualification))
    }
  }

  #Form neat output-----------------------------------------------------------
  origin_allresult<-get_origin_allresult(allresult=allresult,
                                         mingdan=mingdan)
  nbig<-detect_possible_num_cluster(allresult=origin_allresult,
                                    totalnumber=nrow(transformed_datatable),
                                    minrt=minrt)
  if(length(nbig)>1){
    if(bigmark<=1){
      print("The setting of minrt may be too small! The result is meaningless.")
    }
  }
  #Output the results----------------------------------------------------------
  endresult<-list()
  endresult[["origin_allresult"]]<-origin_allresult
  endresult[["number_of_clusters"]]<-nbig
  endresult[["inputparameter"]]<-list()
  endresult[["inputparameter"]][["minrt"]]<-minrt
  endresult[["inputparameter"]][["transformtype"]]<-transformtype
  endresult[["inputparameter"]][["basecluster"]]<-basecluster
  endresult[["inputparameter"]][["emphasizedistance"]]<-emphasizedistance
  endresult[["inputparameter"]][["stoprate"]]<-stop_rate
  endresult[["inputparameter"]][["fanshu"]]<-bn
  endresult[["midparameter"]]<-list()
  endresult[["midparameter"]][["C_jump"]] <- C_jump
  endresult[["midparameter"]][["C_distribution"]] <- C_distribution
  endresult[["midparameter"]][["mediandm"]]<-mediandm
  endresult[["midparameter"]][["C"]]<-C
  endresult[["midparameter"]][["allresult"]]<-allresult
  endresult[["midparameter"]][["densevector"]]<-densevector
  endresult[["midparameter"]][["transitionMatrixOnFirstLevel"]]<-transitionMatrixOnFirstLevel
  stopCluster(cl)
  return(endresult)
}
