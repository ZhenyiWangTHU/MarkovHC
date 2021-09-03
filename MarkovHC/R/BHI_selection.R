#' Recommend basins by biological homogeneity index
#'
#' Function \code{BHI_selection} The function recommends basins
#' by biological homogeneity index.
#' @param SeuratObject The seurat object.
#' @param MarkovObject The MarkovHC object.
#' @param levels An integer vector indicates the customized levels.
#' @param prune A Bloolean parameter indicates wether to merge small basins into
#' big basins. Default is TRUE.
#' @param weed An integer value defines small basins. When the prune parameter is
#' set to TRUE, the basins include less points than weed will be merged into
#' big basins on each level.
#' @param OrgDb OrgDb, and this parameter is passed to 'enrichGO' function in 'clusterProfiler' package.
#' @param ont One of "MF", "BP", and "CC" subontologies. And this parameter is passed to 'enrichGO' function in 'clusterProfiler' package.
#' @param keyType keytype of input gene, and this parameter is passed to 'enrichGO' function in 'clusterProfiler' package.
#' @details This function calculate the biological homogeneity index(BHI) of basins on the customized levels. The definition of
#' BHI is refered to 'clValid' package. And the BHI is in the range [0, 1], with larger values corresponding to more biologically homogeneous basins.
#'
#' To calculate the BHI of a basin, the differentially expressed genes among its sub-basins are found and the BHI is calculated according
#' to the denifition in 'clValid' package. Details are in 'MarkovHC' paper.
#' @return A BHI vector with the names of basins is returned. This function creates a 'GOEachBasin' folder in the current directory and saves the
#' gene ontology enrichment result of each basin in it.
#' @author Zhenyi Wang wangzy17@mails.tsinghua.edu.cn
#' @export

BHI_selection = function(SeuratObject=NULL,
                         MarkovObject=NULL,
                         levels=NULL,
                         prune=TRUE,
                         weed=10,
                         OrgDb='org.Mm.eg.db',
                         ont='all',
                         keyType="ENSEMBL"){
  dir.create('./GOEachBasin')
  # the label matrix
  labels <- fetchLabels(MarkovObject = MarkovObject,
                        MarkovLevels = 1:length(MarkovObject$hierarchicalStructure),
                        prune = FALSE,
                        weed = 10)
  for (i in 2:ncol(labels)) {
    max.temp <- max(table(labels[, i]))
    if (max.temp < weed) {
      next
    }
    else {
      startLv <- i
      break
    }
  }
  labels_temp <- fetchLabels(MarkovObject = MarkovObject,
                             MarkovLevels = startLv:length(MarkovObject$hierarchicalStructure),
                             prune = prune,
                             weed = weed)
  labels[, startLv:ncol(labels)] <- labels_temp[, 1:ncol(labels_temp)]

  for(i in 1:ncol(labels)){
    labels[,i] <- gsub("\\+", "plus", labels[,i])
  }

  BHI_values <- c()
  # enrich the DEGs
  for(i in rev(levels)){
    basins_lvi <- unique(labels[,i])
    if(i==1){break}
    Idents(SeuratObject) <- labels[,(i-1)]
    for(j in basins_lvi){
      basins_lvi_1_in_basinj <- unique( labels[which(labels[,i]==j),(i-1)] )
      if(length(basins_lvi_1_in_basinj)==1){next}
      # find DEGs
      SeuratObject.temp <- subset(SeuratObject, idents = basins_lvi_1_in_basinj)
      SeuratObject.temp.markers <- FindAllMarkers(SeuratObject.temp, only.pos = TRUE, min.pct = 0.1, logfc.threshold = 0.25)
      SeuratObject.temp.markers <- SeuratObject.temp.markers %>% group_by(cluster)
      unique_clusters <- unique(SeuratObject.temp.markers$cluster)
      GO_results <- list()
      for(clusters in unique_clusters){
        DEG_cluster <- subset(SeuratObject.temp.markers, SeuratObject.temp.markers$cluster == clusters)
        if(nrow(DEG_cluster)==0){next}
        # enrich GO
        enrichGO_results <- enrichGO(gene = DEG_cluster$gene,
                                     keyType = keyType,
                                     OrgDb = OrgDb,
                                     ont = ont,
                                     pAdjustMethod = "fdr",
                                     pvalueCutoff = 0.05,
                                     qvalueCutoff  = 0.2,
                                     minGSSize = 3,
                                     maxGSSize = 500,
                                     readable = FALSE)
        enrichGO_results <- as.data.frame(enrichGO_results@result)
        if(nrow(enrichGO_results)==0){next}
        enrichGO_results <- enrichGO_results[order(enrichGO_results[,9], decreasing = TRUE),]
        write.csv(enrichGO_results, file=paste('./GOEachBasin/',
                                               'lv',as.character(i-1),
                                               '_basin',
                                               as.character(clusters),
                                               '.csv',sep=''))
        GO_results <- c(GO_results, list(enrichGO_results))
      }
      GO_terms <- c()
      GO_dataframe <- data.frame()
      for(GOi in 1:length(GO_results)){
        GO_terms <- c(GO_terms, GO_results[[GOi]]$ID)
        GO_dataframe <- rbind(GO_dataframe, GO_results[[GOi]])
      }
      GO_terms <- unique(GO_terms)
      # annotation matrix
      annotation.matrix <- matrix(FALSE, ncol=length(GO_terms), nrow=length(unique(SeuratObject.temp.markers$gene)))
      colnames(annotation.matrix) <- GO_terms
      rownames(annotation.matrix) <- unique(SeuratObject.temp.markers$gene)

      for ( annotation_i in 1:nrow(GO_dataframe) ){
        annot <- GO_dataframe$ID[annotation_i]
        genes <- str_split(GO_dataframe$geneID[annotation_i], pattern = '/')[[1]]
        annotation.matrix[genes, annot] <- TRUE
      }
      # calculate BHI
      statClust <- mapvalues(rownames(annotation.matrix), from=SeuratObject.temp.markers$gene, to=SeuratObject.temp.markers$cluster)
      BHI_basinj <- BHI(statClust=statClust, annotation=annotation.matrix)
      BHI_values <- c(BHI_values, BHI_basinj)
      names(BHI_values)[length(BHI_values)] <- paste('lv',as.character(i),'_basin',as.character(j),sep = '')
    }
  }
  return(BHI_values)
}

BHI <- function(statClust=NULL,
                annotation=NULL) {
  ## initialize BHI vector to 0s
  bhi <- numeric(length(unique(statClust)))
  names(bhi) <- unique(statClust)

  ## for each statClust
  for ( k in unique(statClust) ){
    Ck.bhi <- 0
    Ck.idx <- which(statClust==k) # row indices of this statClust Ck

    if ( length(Ck.idx)<2 ){ next }# only one gene, skip

    ## for each gene in this statClust
    for ( i in Ck.idx ){
      ## ... count how many other genes j in Ck share any (1 or more)
      ## of gene i's annotations:

      B <- which(annotation[i,]==TRUE) # get indices of i's annotations
      if ( length(B)==0 ){ next }# gene i has no annotation

      ## gene's annotations of all other genes j in statClust Ck
      annot <- annotation[Ck.idx[ Ck.idx!= i ],B]
      ## ... add number of genes with at least one shared annotation
      if ( length(B)==1 )      Ck.bhi <- Ck.bhi + sum(annot)
      else if ( length(B) >1 ) Ck.bhi <- Ck.bhi + sum(rowSums(annot)>0)
    }

    nk <- sum(rowSums(annotation[Ck.idx,])>0) # nr. of annot. feat. in Ck
    if ( nk>1 ) bhi[k] <- Ck.bhi / (nk*(nk-1))
  }
  return(mean(bhi, na.rm=TRUE))
}
