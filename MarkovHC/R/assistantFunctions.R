#We deposit assiatant functions for main function "MarkovHC" in this file.
#And would not be exported to users.

##Function1:Transform----------------------------------------------------------
arsinh=function(x=NULL){
 return(log(x+sqrt(x^2+1)))
}

##Function2:Generate the transition probability matrix-------------------------
transition_probability=function(matrix=NULL,
                                densevector=NULL,
                                weightDist=NULL,
                                weightDens=NULL){
 dense<-densevector
 m<-length(densevector)
 G_matrix <- matrix
 for (i in 1:nrow(matrix)) {
   for (j in 1:ncol(matrix)) {
     G_matrix[i,j] <- ((matrix[i,j]^(-weightDist))*((dense[i]/dense[j])^weightDens))
   }
 }
 P <- matrix(data = 0,nrow = nrow(G_matrix), ncol = ncol(G_matrix))
 for (i in 1:nrow(matrix)) {
   min_G_matrix_rowi <- min(G_matrix[i,])
   min_num <- sum(G_matrix[i,]==min_G_matrix_rowi)
   min_index <- which(G_matrix[i,]==min_G_matrix_rowi)
   P[i,min_index] <- 1/min_num
 }
 return(P)
}

##Function3:Calculate the pseudo energy matrix---------------------------------
Calculate_C_Matrix = function(matrix=NULL,
                              densevector=NULL,
                              emphasizedistance=NULL,
                              weightDist=NULL,
                              weightDens=NULL){
 m<-length(densevector)
 dense<-as.numeric(densevector)
 G_matrix <- matrix^emphasizedistance
 for (i in 1:nrow(matrix)) {
   for (j in 1:ncol(matrix)) {
     G_matrix[i,j] <- ((matrix[i,j]^(-weightDist))*((dense[i]/dense[j])^weightDens))
   }
 }
 C <- G_matrix
 for (i in 1:nrow(matrix)) {
   min_G_matrix_rowi <- min(G_matrix[i,])
   for (j in 1:ncol(matrix)) {
     C[i,j] <- G_matrix[i,j]-min_G_matrix_rowi
   }
 }
 return(C)
}

##Function4:Update the transition matrix
update_P = function(C_matrix_updated=NULL,
                    C_cut=NULL){
  # remove zero We can prove it in the next version!
  # cutpoint <- quantile(unique(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)]),probs = C_cut)
  cutpoint <- quantile(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)],probs = C_cut)
  p_updated <- C_matrix_updated
  p_updated[which(C_matrix_updated>cutpoint)] <- 0
  p_updated[which(C_matrix_updated<cutpoint)] <- 1
  p_updated <- p_updated/rowSums(p_updated)
  p_updated_indice <- p_updated
  diag(p_updated_indice) <- 0
  C_cut_step <- C_cut
  while((cutpoint==0)|(all(p_updated_indice==0))){
    C_cut <- C_cut + C_cut_step
    if((C_cut-1)>0){C_cut <- 1}
    #  remove zero We can prove it in the next version!
    #  cutpoint <- quantile(unique(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)]),probs = C_cut)
    cutpoint <- quantile(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)],probs = C_cut)
    p_updated <- C_matrix_updated
    p_updated[which(C_matrix_updated>cutpoint)] <- 0
    p_updated[which(C_matrix_updated<cutpoint)] <- 1
    p_updated <- p_updated/rowSums(p_updated)
    p_updated_indice <- p_updated
    diag(p_updated_indice) <- 0
  }
  #print(C_matrix_updated)
  #print(p_updated)
  #print(paste('cutpoint is ', as.character(cutpoint),'.', sep = ''))
  return(list(p_updated, cutpoint))
}

#Function5:Arrange the coordinate of the hierarchical structure
arr_coordinate = function(MarkovObject=NULL,
                          startLevel=NULL,
                          endLevel=NULL){
  orderedMarkov <- orderMarkovHC(MarkovObject=MarkovObject,
                                 MarkovLevels=1:length(MarkovObject[["hierarchicalStructure"]]),
                                 orderLevels=1:length(MarkovObject[["hierarchicalStructure"]]))

  firstOrder <- unique(orderedMarkov[,startLevel+1])
  firstCoordinate <- vector(mode = 'integer',length = length(firstOrder))
  index_temp <- 1
  for (i in firstOrder) {
    firstCoordinate[i] <- index_temp
    index_temp <- index_temp+1
  }
  coordinate <- list()
  coordinate[[1]] <- firstCoordinate
  coordinateIndex <- 2
  for (i in (startLevel+1):endLevel) {
    coordinate_temp <- vector(mode = 'integer',
                              length = length(MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(i),sep='')]][['basins']]))
    basinSize <- MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(i-1),sep='')]][['basinPoints']]%>%lengths()

    for (j in 1:length(coordinate_temp)) {
      basinsIntemp <- MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(i),sep='')]][['basins']][[j]]
      biggestBasin <- basinsIntemp[which.max(basinSize[basinsIntemp])[1]]
      coordinate_temp[j] <- coordinate[[coordinateIndex-1]][biggestBasin]
    }

    coordinate[[coordinateIndex]] <- coordinate_temp
    coordinateIndex <- coordinateIndex+1
  }

  for (i in 1:length(coordinate)) {
    coordinate[[i]] <- cbind(coordinate[[i]],rep((1+0.5*(i-1)),length(coordinate[[i]])))
  }

  return(coordinate)
}
