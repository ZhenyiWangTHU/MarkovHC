#' @export
#' plot the hierarchical structure
#' The color vector customize the colors of basins on the customized lowest level.
plotHierarchicalStructure = function(MarkovObject=NULL,
                                     MarkovLevels=NULL,
                                     colorVector=NULL){
    #Build a tree
    options(stringsAsFactors = TRUE)
    clusterings <- fetchLabels(MarkovObject,MarkovLevels)
    for(i in 1:ncol(clusterings)){
      for(j in 1:nrow(clusterings)){
        clusterings[j,i] <- str_split(clusterings[j,i],'\\+')[[1]][1]
        }
    }
    nodes <- MarkovHC::getNodes(clusterings)
    #order the nodes
    nodes <- as.data.frame(nodes)
    nodes$level <- as.numeric(nodes$level)
    nodes$basin <- as.numeric(nodes$basin)
    nodes <- doBy::orderBy(as.formula(~level+basin), nodes)
    rownames(nodes) <- 1:nrow(nodes)
    #需要指定nodes的顺序，否则颜色会乱掉
    nodes$Node <- factor(nodes$Node, levels = nodes$Node)
    treeEdge <- fetchLabels(MarkovObject,MarkovLevels)
    for(i in 1:ncol(treeEdge)){
      for(j in 1:nrow(treeEdge)){
        treeEdge[j,i] <- str_split(treeEdge[j,i],'\\+')[[1]][1]
      }
    }
    colnames(treeEdge) <- paste('level',MarkovLevels,sep='')
    edges <- MarkovHC::getEdges(treeEdge)

    #graph
    graph <- edges %>%
    # Remove edges without any cell...
    filter(TransCount > 0) %>%
    # Rename the nodes
    mutate(FromNode = paste0("level", FromRes, "_", FromClust)) %>%
    mutate(ToNode = paste0("level", ToRes, "_", ToClust)) %>%
    # Reorder columns
    dplyr::select(FromNode, ToNode, everything()) %>%
    # Build a graph using igraph
    graph_from_data_frame(vertices = nodes)
    #arrange the coordinate
    nodesPosition <- arr_coordinate(MarkovObject = MarkovObject,
                                    startLevel = MarkovLevels[1],
                                    endLevel = MarkovLevels[length(MarkovLevels)])
    #print(nodesPosition)
    nodesPosition <- do.call(rbind,lapply(nodesPosition, data.frame))
    #print(nodesPosition)
    #arrange the colors
    arrangedColorVector <- nodesPosition
    arrangedColorVector$color <- character(length = nrow(nodesPosition))

    #colorVector <- colorVector[nodesPosition[1:length(MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(MarkovLevels[1]),sep='')]][["basins"]]),1]]
    #print(nodesPosition[1:length(MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(MarkovLevels[1]),sep='')]][["basins"]]),1])
    #print(colorVector)
    position_lowestLevel <- nodesPosition[1:length(MarkovObject[["hierarchicalStructure"]][[paste('level',as.character(MarkovLevels[1]),sep='')]][["basins"]]),1]
    for (i in 1:length(position_lowestLevel)) {
      arrangedColorVector[which(arrangedColorVector[,1]==position_lowestLevel[i]),3] <- colorVector[i]
    }
    #print(arrangedColorVector)
    arrangedColorVector <- arrangedColorVector$color
    colnames(nodesPosition) <- c('x','y')

    #for [Package ggraph version 2.0.2 ]
    #ggraph(graph,
    #       layout = "manual",
    #       x=nodesPosition[,1],
    #       y=nodesPosition[,2],
    #       circular = FALSE)
    ggraph(graph,
           layout = "manual",
           node.positions=nodesPosition,
           circular = FALSE) +
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

}
