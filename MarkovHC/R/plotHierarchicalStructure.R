#' Plot the MarkovHC hierarchical structure or return the tree
#'
#' Function \code{plotHierarchicalStructure} The function plots the MarkovHC hierarchical structure or returns the tree.
#' @param MarkovObject The output of the function, \code{MarkovHC}.
#' @param MarkovLevels An integer vector indicates the customized levels.
#' @param colorVector A character vector set the colors of basins on the first customized level.
#' The length of this vector must equal to the number of basins on the first customized level.
#' @param plot A Bloolean parameter indicates wether to plot the MarkovHC hierarchical structure.
#' If this parameter was set to FALSE, this function will not plot and a phylo tree will be returned. The
#' phylo tree could be plotted by some tree visualization packages, such as "ggtree".
#' @param prune A Bloolean parameter indicates wether to merge small basins into
#' big basins. Default is TRUE.
#' @param weed An integer value defines small basins. When the prune parameter is
#' set to TRUE, the basins include less points than weed will be merged into
#' big basins on each level. Default is 10.
#' @details This function plots the MarkovHC hierarchical structure or returns the tree.
#' @return When the plot parameter is set to TRUE, this function plots the MarkovHC hierarchical structure. When the plot
#' parameter is set to FALSE, this function returns a phylo tree and it could be plotted by other tree visualization
#' packages, such as "ggtree".
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export
plotHierarchicalStructure = function(MarkovObject=NULL,
                                     MarkovLevels=NULL,
                                     colorVector=NULL,
                                     plot=TRUE,
                                     prune=TRUE,
                                     weed=10){
  #Build a tree
  options(stringsAsFactors = TRUE)
  if(prune==FALSE){weed <- 0}
  clusterings <- fetchLabels(MarkovObject,MarkovLevels, prune=prune, weed=weed)
  treeEdge <- clusterings
  nodes <- getNodes(clusterings)

  #order the nodes
  nodes <- as.data.frame(nodes)
  nodes$level <- as.numeric(nodes$level)
  nodes$basin <- as.factor(nodes$basin)
  nodes <- doBy::orderBy(as.formula(~level+basin), nodes)
  rownames(nodes) <- 1:nrow(nodes)

  #需要重排nodes的顺序，否则坐标和颜色会乱掉
  orderLevels_formula <- MarkovLevels[order(MarkovLevels, decreasing = TRUE)]%>%
    paste('lv',., sep='')%>%
    paste(., collapse = '+')%>%
    paste('~',., sep = '')
  orderedMarkov <- doBy::orderBy(as.formula(orderLevels_formula), clusterings)
  nodesOrder <- c()
  for(i in 1:ncol(orderedMarkov)){
    nodesOrder <- c(nodesOrder, paste(colnames(orderedMarkov)[i],unique(orderedMarkov[,i]),sep='_'))
  }

  nodes$Node <- factor(nodes$Node, levels = nodesOrder)
  nodes <- nodes[order(nodes$Node, decreasing = FALSE),]

  edges <- getEdges(treeEdge)

  if(plot==FALSE){
    edges$edge.length <- edges$FromRes
    for(i in 2:length(MarkovLevels)){
      index_temp <- which((edges$FromRes==MarkovLevels[i-1])&(edges$ToRes==MarkovLevels[i]))
      edges[index_temp,"edge.length"] <- MarkovObject$hierarchicalStructure[[MarkovLevels[i-1]]]$energyCutpoint
    }
  }

  #graph
  edges <- edges %>%
    # Remove edges without any cell...
    #filter(TransCount > 0) %>%
    dplyr::filter(TransCount > weed) %>%
    # Rename the nodes
    mutate(FromNode = paste0("lv", FromRes, "_", FromClust)) %>%
    mutate(ToNode = paste0("lv", ToRes, "_", ToClust)) %>%
    # Reorder columns
    dplyr::select(FromNode, ToNode, everything())

  graph <- edges%>%
    # Build a graph using igraph
    graph_from_data_frame(vertices = nodes)

  if(plot==TRUE){
    #arrange the coordinate
    if(prune==TRUE){
      clusterings.temp <- fetchLabels(MarkovObject,MarkovLevels, prune=FALSE, weed=weed)
      orderedMarkov.temp <- doBy::orderBy(as.formula(orderLevels_formula), clusterings.temp)
      nodesPosition.temp <- arr_coordinate_labelMatrix(orderedMarkov=orderedMarkov.temp)
      nodesPosition.temp <- do.call(rbind,lapply(nodesPosition.temp, data.frame))
      nodesOrder.temp <- c()
      for(i in 1:ncol(orderedMarkov.temp)){
        nodesOrder.temp <- c(nodesOrder.temp, paste(colnames(orderedMarkov.temp)[i],unique(orderedMarkov.temp[,i]),sep='_'))
      }
      rownames(nodesPosition.temp) <- nodesOrder.temp
      nodesPosition.temp <- subset(nodesPosition.temp, rownames(nodesPosition.temp)%in%nodes$Node)
      nodesPosition.temp$nodes <- factor(rownames(nodesPosition.temp), levels = nodes$Node)
      nodesPosition.temp <- nodesPosition.temp[order(nodesPosition.temp$nodes, decreasing = FALSE),]
      nodesPosition.temp <- nodesPosition.temp[,-which(colnames(nodesPosition.temp)=='nodes')]
      nodesPosition <- nodesPosition.temp
    }else{
      nodesPosition <- arr_coordinate_labelMatrix(orderedMarkov=orderedMarkov)
      nodesPosition <- do.call(rbind,lapply(nodesPosition, data.frame))
    }
    #arrange the colors
    arrangedColorVector <- nodesPosition
    arrangedColorVector$color <- character(length = nrow(nodesPosition))

    position_lowestLevel <- nodesPosition[1:length(unique(clusterings[,1])),1]
    for (i in 1:length(position_lowestLevel)) {
      arrangedColorVector[which(arrangedColorVector[,1]==position_lowestLevel[i]),3] <- colorVector[i]
    }

    arrangedColorVector <- arrangedColorVector$color
    colnames(nodesPosition) <- c('x','y')

    #for [Package ggraph version 2.0.2 ]
    ggraph(graph,
           layout = "manual",
           x=nodesPosition[,1],
           y=nodesPosition[,2],
           circular = FALSE)+
      #for [Package ggraph version < 2.0.2 ]
      #     ggraph(graph,
      #            layout = "manual",
      #            node.positions=nodesPosition,
      #            circular = FALSE) +
      # Plot the edges, colour is the number of cells and transparency is the
      # proportion contribution to the new cluster
      geom_edge_link(arrow = arrow(length = unit(7, 'mm')),
                     end_cap = circle(10, "mm"),
                     edge_width = 2
                     #aes(colour = log(TransCount), alpha = TransPropTo)
      ) +
      # Plot the nodes, size is the number of cells
      geom_node_point(aes(colour = factor(nodes$Node),
                          size = Size)) +
      scale_color_manual(
        values = arrangedColorVector)+
      geom_node_text(aes(label = basin), size = 10,
                     family = "sans",
                     color = "black",
                     fontface = "bold") +
      xlim(min(nodesPosition[,1])-0.2,max(nodesPosition[,1])+0.2) +
      ylim(min(nodesPosition[,2])-0.2,max(nodesPosition[,2])+0.2) +
      # Adjust the scales
      scale_size(range = c(10, 35)) +
      #scale_edge_colour_gradientn(colours = viridis(100)) +
      # Add legend labels
      guides(colour = guide_legend(title = "basin colors",
                                   title.position = "top", order = 1),
             size = guide_legend(title = "basin size", title.position = "top", order = 2)
             # edge_alpha = guide_legend(title = "basin prop",
             #                           title.position = "top", nrow = 2, order = 4),
             # edge_colour = guide_edge_colorbar(title = "cell count (log)",
             #                                   title.position = "top", order = 3)
      ) +
      # Remove the axes as they don't really mean anything
      theme_void() +
      theme(legend.position = "right",
            legend.text = element_text(size=20,
                                       family = "sans",
                                       color = "black",
                                       face = "bold"),
            legend.title = element_text(size=25,
                                        family = "sans",
                                        color = "black",
                                        face = "bold")
      )
  }else{
    print('Return a tree for plot.')
    tree <- treeio::as.phylo(graph)
    tree <- treeio::as_tibble(tree)
    tree$branch.length <- as.numeric(mapvalues(tree$label, from=edges$FromNode,to=edges$edge.length))
    tree <- treeio::as.phylo(tree)
    return(tree)
  }
}
