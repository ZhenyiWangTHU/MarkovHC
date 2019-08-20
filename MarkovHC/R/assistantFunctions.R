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
  cutpoint <- quantile(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)],probs = C_cut)
  p_updated <- C_matrix_updated
  p_updated[which(C_matrix_updated>cutpoint)] <- 0
  p_updated[which(C_matrix_updated<cutpoint)] <- 1
  p_updated <- p_updated/rowSums(p_updated)
  p_updated_indice <- p_updated
  diag(p_updated_indice) <- 0
  while((cutpoint==0)|(all(p_updated_indice==0))){
    C_cut <- C_cut + 0.01
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





