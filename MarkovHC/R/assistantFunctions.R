# We deposit assiatant functions in this file.
# And would not be exported to users.

# Calculate the pseudo energy matrix---------------------------------
transition_probability_pseudo_energy_matrix=function(matrix=NULL, D.matrix=NULL){
  matrix <- as.matrix(matrix)
  G_matrix <- (matrix^(-2))*D.matrix
  P <- matrix(data = 0,nrow = nrow(G_matrix), ncol = ncol(G_matrix))
  C <- matrix(data = 0,nrow = nrow(G_matrix), ncol = ncol(G_matrix))
  for (i in 1:nrow(matrix)) {
    min_G_matrix_rowi <- min(G_matrix[i,])
    min_num <- sum((abs(G_matrix[i,]-min_G_matrix_rowi)<exp(-30)))
    min_index <- which((abs(G_matrix[i,]-min_G_matrix_rowi)<exp(-30))==TRUE)
    P[i, min_index] <- 1/min_num
    C[i, ] <- G_matrix[i,]-min_G_matrix_rowi
  }
  return(list(P=P,C=C))
}

# Get a vector to judge recurrent state-----------------------------
judge_RS<-function(P=NULL){
  RSmatrix<-eigen(t(P),symmetric = FALSE)
  #the eigen value is close to 1 enough
  eigen_value_one<-which(abs(RSmatrix$values-1)<exp(-30))
  if(length(eigen_value_one)==1){
    return(abs(RSmatrix$vectors[,eigen_value_one]))
  }else{
    RS<-apply(abs(RSmatrix$vectors[,eigen_value_one]),1,sum)
  }
  return(RS)
}

# Update the transition matrix--------------------------------------
update_P = function(C_matrix_updated=NULL,
                    C_cut=NULL){
  cutpoint <- quantile(unique(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)]), probs = C_cut)
  p_updated <- C_matrix_updated
  p_updated[which(C_matrix_updated>cutpoint)] <- 0
  p_updated[which(C_matrix_updated<=cutpoint)] <- 1
  diag(p_updated) <- 1
  p_updated <- p_updated/rowSums(p_updated)
  p_updated_indice <- p_updated
  diag(p_updated_indice) <- 0
  C_cut_step <- C_cut
  while(all(p_updated_indice==0)){
    C_cut <- C_cut + C_cut_step
    #overcome R accuracy
    if(((C_cut-1)>0)&((C_cut-1)<C_cut_step)){C_cut <- 1}else if((C_cut-1)>=C_cut_step){break}
    cutpoint <- quantile(unique(C_matrix_updated[which(is.infinite(C_matrix_updated)==FALSE)]), probs = C_cut)
    p_updated <- C_matrix_updated
    p_updated[which(C_matrix_updated>cutpoint)] <- 0
    p_updated[which(C_matrix_updated<=cutpoint)] <- 1
    diag(p_updated) <- 1
    p_updated <- p_updated/rowSums(p_updated)
    p_updated_indice <- p_updated
    diag(p_updated_indice) <- 0
  }
  return(list(p_updated, cutpoint))
}

# get Edges---------------------------------------------------------
getEdges <- function(clusterings=NULL) {

  # Loop over the different resolutions
  transitions <- lapply(1:(ncol(clusterings) - 1), function(i) {

    # Extract two neighbouring clusterings
    from.res <- colnames(clusterings)[i]
    to.res <- colnames(clusterings)[i + 1]

    # Get the cluster names
    from.clusters <- sort(unique(clusterings[, from.res]))
    to.clusters <- sort(unique(clusterings[, to.res]))

    # Get all possible combinations
    trans.df <- expand.grid(FromClust = from.clusters,
                            ToClust = to.clusters,
                            stringsAsFactors = FALSE)

    # Loop over the possible transitions
    trans <- apply(trans.df, 1, function(x) {
      from.clust <- x[1]
      to.clust <- x[2]

      # Find the cells from those clusters
      is.from <- clusterings[, from.res] == from.clust
      is.to <- clusterings[, to.res] == to.clust

      # Count them up
      trans.count <- sum(is.from & is.to)

      # Get the sizes of the two clusters
      from.size <- sum(is.from)
      to.size <- sum(is.to)

      # Get the proportions of cells moving along this edge
      trans.prop.from <- trans.count / from.size
      trans.prop.to <- trans.count / to.size

      return(c(trans.count, trans.prop.from, trans.prop.to))
    })

    # Tidy up the results
    trans.df$FromRes <- as.numeric(gsub("lv", "", from.res))
    trans.df$ToRes <- as.numeric(gsub("lv", "", to.res))
    trans.df$TransCount <- trans[1, ]
    trans.df$TransPropFrom <- trans[2, ]
    trans.df$TransPropTo <- trans[3, ]

    return(trans.df)
  })

  # Bind the results from the different resolutions together
  transitions <- do.call("rbind", transitions)

  # Tidy everything up
  #levs <- sort(as.numeric(levels(transitions$ToClust)))

  #levs <- sort(levels(transitions$FromClust))
  #transitions <- transitions %>%
  #  mutate(FromClust = factor(FromClust,
  #                            levels = levs))  %>%
  #  mutate(ToClust = factor(ToClust, levels = levs))

  return(transitions)
}

# get Nodes--------------------------------------------------------------------
getNodes <- function(clusterings=NULL) {
  nodes <- clusterings %>%
    tidyr::gather(key = level, value = basin) %>%
    group_by(level, basin) %>%
    dplyr::summarise(Size = n()) %>%
    dplyr::ungroup() %>%
    mutate(level = stringr::str_replace(level, "lv", "")) %>%
    mutate(level = level, basin = basin) %>%
    mutate(Node = paste0("lv", level ,'_', basin))%>%
    dplyr::select(Node,everything())
}

# Arrange the coordinate of the hierarchical structure-------------------------
arr_coordinate_labelMatrix = function(orderedMarkov=NULL){
  orderedMarkov <- orderedMarkov
  firstOrder <- unique(orderedMarkov[,1])
  firstCoordinate <- 1: length(firstOrder)

  coordinate <- list()
  coordinate[[1]] <- firstCoordinate
  coordinateIndex <- 2
  for (i in 2:ncol(orderedMarkov)) {
    coordinate_temp_i <- vector(mode = 'integer',
                                length = length(unique(orderedMarkov[,i])))
    basin_temp_i <- unique(orderedMarkov[,i])
    basin_temp_i_1 <- unique(orderedMarkov[,(i-1)])

    for (j in 1:length(coordinate_temp_i)) {
      pointsIntemp_i_1 <- orderedMarkov[which(orderedMarkov[,i]==basin_temp_i[j]), (i-1)]

      maxBasin <- names(table(pointsIntemp_i_1))[table(pointsIntemp_i_1) == max(table(pointsIntemp_i_1))]

      maxBasin_index <- which(basin_temp_i_1==maxBasin)

      coordinate_temp_i[j] <- coordinate[[coordinateIndex-1]][maxBasin_index]
    }
    coordinateIndex <- coordinateIndex+1
    coordinate <- c(coordinate, list(coordinate_temp_i))
  }

  for (i in 1:length(coordinate)) {
    coordinate[[i]] <- cbind(coordinate[[i]],rep((1+0.5*(i-1)),length(coordinate[[i]])))
  }

  return(coordinate)
}

# find peaks--------------------------------------------------------
# refer to https://github.com/stas-g/findPeaks
find_peaks = function (x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  pks
}
