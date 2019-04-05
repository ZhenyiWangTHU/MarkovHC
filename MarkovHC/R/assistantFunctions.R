#We deposit assiatant functions for main function "MarkovHC" in this file.
#And would not export them to users.
#Transform---------------------------------------------------------------------
arsinh<-function(x=NULL){
 return(log(x+sqrt(x^2+1)))
}

#Calculate the density vector--------------------------------------------------
get_densityvector=function(matrix=NULL){
 distancematrix <- dist(matrix, method = "euclidean", diag = F, upper = T, p = 2)
 distancematrix <- as.matrix(distancematrix)
 dv<-as.vector(distancematrix)
 radius <- quantile(dv,0.05)
 densityvector <- numeric(nrow(distancematrix))
 for (i in 1:nrow(distancematrix)) {
   densityvector[i] <- sum(distancematrix[i,]<=radius)
 }
 return(densityvector)
}

#Generate a transition matrix on the first level------------------------------
# transition_first_level=function(matrix=NULL,densevector=NULL,lambda=NULL){
#   dense<-densevector
#   m<-length(densevector)
#   pinfpar<-function(i,matrix_temp=matrix,dense_temp=dense,lambda_temp=lambda,m_temp=m){
#     transi<-(matrix_temp[i,]==0)
#     #here is a little different from paper,sqrt the density ratio
#     ddense<-lambda_temp*sqrt(dense_temp/dense_temp[i])
#     m1<-(ddense==1)
#     m2<-(ddense>1)
#     a1<-which((m1*transi)==1)
#     a2<-which((m2*transi)==1)
#     p<-rep(0,m_temp)
#     p[a1]<-(1-exp(-1))*sqrt(dense_temp[a1]/dense_temp[i])
#     p[a2]<-sqrt(dense_temp[a2]/dense_temp[i])
#     p<-p/sum(p)
#   }
#   Pinf<-foreach(x=1:m,.combine='rbind') %dopar% pinfpar(x)
#   return(Pinf)
# }
transition_first_level=function(matrix=NULL,densevector=NULL,weightDist=NULL,weightDens=NULL){
 dense<-densevector
 m<-length(densevector)
 G_matrix <- matrix
 for (i in 1:nrow(matrix)) {
   for (j in 1:ncol(matrix)) {
     G_matrix[i,j] <- ((matrix[i,j]^weightDist)*(dense[i]/dense[j])^weightDens)
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

#Generate pseudo potential difference matrix C--------------------------------
# Calculate_C_Matrix<-function(GD=NULL,densevector=NULL,lambda=NULL, eta=NULL,emphasizedistance=NULL){
#   m<-length(densevector)
#   dense<-as.numeric(densevector)
#   C<-(GD^emphasizedistance)*eta
#   C<-C*sqrt(dense)
#   C<-t(t(C)*sqrt((1/dense)))
#   for(i in 1:m){
#    transi<-which(lambda*sqrt(dense/dense[i])<1)
#    C[i,transi]<-C[i,transi]+log(sqrt(dense[i]/dense[transi])/lambda)
#   }
#   return(C)
# }
Calculate_C_Matrix<-function(GD=NULL,densevector=NULL, emphasizedistance=NULL,weightDist=NULL,weightDens=NULL){
 m<-length(densevector)
 dense<-as.numeric(densevector)
 G_matrix <- GD^emphasizedistance
 for (i in 1:nrow(GD)) {
   for (j in 1:ncol(GD)) {
     G_matrix[i,j] <- ((GD[i,j])^weightDist)*((dense[i]/dense[j])^weightDens)
   }
 }
 C <- G_matrix
 for (i in 1:nrow(GD)) {
   min_G_matrix_rowi <- min(G_matrix[i,])
   min_num <- sum(G_matrix[i,]==min_G_matrix_rowi)
   for (j in 1:ncol(GD)) {
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










