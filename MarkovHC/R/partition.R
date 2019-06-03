#' Partition the state space
#' 
partition = function(C_matrix=NULL, RS_vector=NULL){
  processed_indice <- integer(length = length(RS_vector))
  recurrent_indice <- which(RS_vector>0)
  C_graph <- graph_from_adjacency_matrix(adjmatrix=C_matrix, mode="directed", weighted=TRUE, diag=TRUE)
  
  for(i in recurrent_indice){
    if(processed_indice[i]==1){next}
    
  }
  
}