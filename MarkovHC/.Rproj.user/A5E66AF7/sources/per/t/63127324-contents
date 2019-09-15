#' @export
#get Edges
getEdges <- function(clusterings=NULL) {

  # Loop over the different resolutions
  transitions <- lapply(1:(ncol(clusterings) - 1), function(i) {

    # Extract two neighbouring clusterings
    from.res <- sort(colnames(clusterings))[i]
    to.res <- sort(colnames(clusterings))[i + 1]

    # Get the cluster names
    from.clusters <- sort(unique(clusterings[, from.res]))
    to.clusters <- sort(unique(clusterings[, to.res]))

    # Get all possible combinations
    trans.df <- expand.grid(FromClust = from.clusters,
                            ToClust = to.clusters)

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
    trans.df$FromRes <- as.numeric(gsub("level", "", from.res))
    trans.df$ToRes <- as.numeric(gsub("level", "", to.res))
    trans.df$TransCount <- trans[1, ]
    trans.df$TransPropFrom <- trans[2, ]
    trans.df$TransPropTo <- trans[3, ]

    return(trans.df)
  })

  # Bind the results from the different resolutions together
  transitions <- do.call("rbind", transitions)

  # Tidy everything up
  #levs <- sort(as.numeric(levels(transitions$ToClust)))
  levs <- sort(levels(transitions$FromClust))
  transitions <- transitions %>%
    mutate(FromClust = factor(FromClust,
                              levels = levs))  %>%
    mutate(ToClust = factor(ToClust, levels = levs))

  return(transitions)
}

#' @export
# get Nodes
getNodes <- function(clusterings=NULL) {
  nodes <- clusterings %>%
    tidyr::gather(key = level, value = basin) %>%
    group_by(level, basin) %>%
    summarise(Size = n()) %>%
    ungroup() %>%
    mutate(level = stringr::str_replace(level, "level", "")) %>%
    mutate(level = level, basin = basin) %>%
    mutate(Node = paste0("level", level ,'_', basin))%>%
    dplyr::select(Node, everything())
}
















