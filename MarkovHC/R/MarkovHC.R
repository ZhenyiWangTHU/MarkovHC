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
#' @param basecluster A character parameter indicating what kind of simple
#' clustering methods would be used as a preliminary processinng method.
#' Available choices include \code{single}, \code{complete}, \code{average}
#' and \code{kmeans}. Default is 'average'.
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
              basecluster="average",
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

  #basecluster
  if(!is.character(basecluster)){
    print("The type of 'basecluster' should be character!")
    return(NULL)
  }
  if((basecluster=="single")|(basecluster=="complete")|
     (basecluster=="average")|(basecluster=="kmeans")){
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

  #Do parallel-----------------------------------------------------------------
  ncore<-detectCores()
  cl <- makeCluster(getOption("cl.cores", ncore))
  registerDoParallel(cl)

  ##step02.do preclustering----------------------------------------------------
  if(dobasecluster==TRUE){
    #do clustering on the first level
    #Use one type of hierarchical clustering as the basic clustering tool
    if((basecluster=="single")|(basecluster=="complete")|(basecluster=="average")){
      hresult<-hclust(dist(transformed_matrix,method = "minkowski",p=bn),method = basecluster)
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
        hresult_cut <- cutree(hresult,k=baseclusternum)
      }else{
        #user sets the baseclusternum parameter
        hresult_cut <- cutree(hresult,k=baseclusternum)
      }
    }else{
      if (is.null(baseclusternum)){
        baseclusternum <- ceiling(nrow(transformed_matrix)/10)
        kmeansresult <- kmeans(transformed_matrix, centers=baseclusternum,iter.max = 500)
        hresult_cut <- kmeansresult$cluster
      }else{
        #user sets the baseclusternum parameter
        kmeansresult <- kmeans(transformed_matrix, centers=baseclusternum,iter.max = 500)
        hresult_cut <- kmeansresult$cluster
      }
    }
  }else{
    #do not do clustering on the first level
    hresult_cut <- 1:nrow(transformed_matrix)
  }

  ##step03.estimate the density of each state----------------------------------
  densevector <-

  ##step04.calculate the distance matrix
  dm <- as.matrix(dist(statematrix,method = "minkowski",p=bn))


  ##step05. find recurrent classes on the first level--------------------------
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
                                    totalnumber=nrow(transformed_matrix),
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