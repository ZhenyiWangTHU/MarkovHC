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
Calculate_C_Matrix<-function(matrix=NULL,
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





#Union of a part of list components-------------------------------------------
bing<-function(c,m){
 ## c is a list containing several sets and m is an index set, we will
 ## union those sets in c whose index is in m
 g<-vector()
 for(i in m){
   g<-union(g,c[[i]])
 }
 return(g)
}


#An assistant function for calculating cluster Ratio--------------------------
clusterRatio = function(endResult = None,
                       level=None,
                       clusterNum=None){
 clusterSize <- lengths(endResult$origin_allresult[[level]][[2]])
 clusterSize <- sort(clusterSize)
 ratio = sum(clusterSize[(length(clusterSize)-clusterNum+1):length(clusterSize)])/sum(clusterSize)
 return(ratio)
}

#An assistant function to modify the first transition matrix
modify_P<-function(P_tran,record_cut,densevector){
 #this function is used to strengthen the status of center points
 if(length(record_cut)==1){
   return(P_tran)
 }
 P<-P_tran
 for(i in 1:(length(record_cut)-1)){
   P[(record_cut[i]+1):(record_cut[i+1]-1),record_cut[i]]<-1
   startp<-(record_cut[i]+1)
   endp<-(record_cut[i+1]-1)
   s<-which(densevector[startp:endp]/densevector[record_cut[i]]>1/1.05)
   P[record_cut[i],(startp:endp)[s]]<-1
 }
 P<-(P/rowSums(P))
 return(P)
}

#An assistant function to modify the C matrix
modify_C<-function(C,record_cut,densevector){
 #this function is used to do some resonable modifications on C to get rid of the influence of taking representatives
 if(length(record_cut)==1){
   return(C)
 }
 CT<-C
 for(i in 1:(length(record_cut)-1)){
   CT[(record_cut[i]+1):(record_cut[i+1]-1),record_cut[i]]<-0
   startp<-(record_cut[i]+1)
   endp<-(record_cut[i+1]-1)
   s<-which(densevector[startp:endp]/densevector[record_cut[i]]>1/1.05)
   CT[record_cut[i],(startp:endp)[s]]<-0
 }
 return(CT)
}










