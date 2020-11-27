# use spaltter to generate simulation data


```R
library(splatter)
library(scater)
library(ggplot2)
library(Seurat)
library(car)
library(rgl)
library(MarkovHC)
options(rgl.useNULL=FALSE)
# Linear paths
params.groups <- newSplatParams(batchCells = 1000, 
                                nGenes = 5000,
                                group.prob = c(0.35, 0.3,0.35),
                                de.prob = 0.5, de.facLoc = 0.2,
                                path.from = c(0, 1, 2),
                                path.skew = c(1,0,0),
                                path.length = c(10,10,10))
sim2Object <- splatSimulatePaths(params.groups,verbose = FALSE)
sim2 <- scater::normalize(sim2Object)
sim2_plot <- scater::plotPCA(sim2, colour_by = "Group") + ggtitle("Linear paths")
sim2_plot
mat2<-counts(sim2Object)
```


```R
#Figure theme
mytheme <- theme(panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(size = 1,
                                          colour = "black"),
                 axis.title.x =element_text(size=20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold"),
                 axis.text.x = element_text(size = 20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold",
                                            vjust = 0,
                                            hjust = 0),
                 axis.text.y = element_text(size = 20,
                                            family = "sans",
                                            color = "black",
                                            face = "bold",
                                            vjust = 0,
                                            hjust = 1),
                 axis.title.y=element_text(size=20,
                                           family = "sans",
                                           color = "black",
                                           face = "bold"),
                 legend.text = element_text(size=15,
                                            family = "sans",
                                            color = "black",
                                            face = "bold"),
                 legend.title = element_text(size=15,
                                             family = "sans",
                                             color = "black",
                                             face = "bold"),
                 legend.background = element_blank(),
                 legend.key = element_blank()
)
```

# use Seurat to preprocess data


```R
PATHobject <- CreateSeuratObject(counts = mat2,
                                 project = 'PATH',
                                 min.cells = 10,
                                 min.feature = 50)
PATHobject@meta.data$Group <- sim2Object$Group
VlnPlot(PATHobject, features = c("nFeature_RNA", "nCount_RNA"), ncol = 2)
PATHobjectplot <- FeatureScatter(PATHobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
PATHobjectplot

PATHobject <- NormalizeData(PATHobject, normalization.method = "LogNormalize", scale.factor = 10000)

PATHobject <- FindVariableFeatures(PATHobject, selection.method = "vst", nfeatures = 3000)
# Identify the 10 most highly variable genes
PATHobjecttop10 <- head(VariableFeatures(PATHobject), 10)
# plot variable features with and without labels
PATHobjectplot1 <- VariableFeaturePlot(PATHobject)
PATHobjectplot2 <- LabelPoints(plot = PATHobjectplot1, points = PATHobjecttop10, repel = TRUE)
PATHobjectplot2
PATHobject <- ScaleData(PATHobject, features = rownames(PATHobject))
PATHobject <- RunPCA(PATHobject, features = VariableFeatures(object = PATHobject), verbose=FALSE)
ElbowPlot(PATHobject, ndims = 50)
```

# run MarkovHC


```R
MarkovHC_PATH <- MarkovHC(origin_matrix=Embeddings(object = PATHobject, reduction = "pca")[,1:5]%>%t(),
                          transformtype="none",
                          KNN=50,
                          basecluster="kmeans",
                          dobasecluster=TRUE,
                          baseclusternum=200,
                          emphasizedistance=1,
                          weightDist=2,
                          weightDens=0.5,
                          cutpoint=0.01,
                          showprocess=FALSE,
                          bn=2,
                          minBasinSize=0.2,
                          noiseBasinSize=20)
```


```R
layout <- sim2_plot$data%>%as.data.frame()
layout$stages <- layout$colour_by
pdf(file = './PATH.A.pdf', width = 9, height = 5.5)
ggplot(data=layout, mapping =  aes(x=X, y=Y)) +
  geom_point(size=2, shape=21, aes(fill=stages), color="#525252")+
  xlim(min(layout$X)-1,max(layout$X)+1)+
  ylim(min(layout$Y)-1,max(layout$Y)+1)+
  mytheme
dev.off()

#Plot the points on PC1,2 and 3 
sim2_PCA <- scater::runPCA(sim2, ncomponents = 5)
sim2_PCA_cor <- sim2_PCA@reducedDims$PCA

pdf('./PATH.PC1.pdf')
g <- ggplot(as.data.frame(sim2_PCA_cor), aes(PC1))
g + geom_density(aes(fill=factor(layout$stages)), alpha=0.5) + 
  #geom_vline(xintercept=c(-6.8,-6.8), linetype="longdash", size=1)+
  labs(title="", 
       x="PC1",
       fill="stage") +
  scale_fill_manual(values = c("Path1" = alpha("#662506",0.2),
                               'Path2'=alpha("#006d2c",0.2),
                               'Path3'= alpha("#08519c",0.2)),
                    breaks = c("Path1","Path2",'Path3'))+
  mytheme
dev.off()

pdf('./PATH.PC2.pdf')
g <- ggplot(as.data.frame(sim2_PCA_cor), aes(PC2))
g + geom_density(aes(fill=factor(layout$stages)), alpha=0.5) + 
  #geom_vline(xintercept=c(-6.8,-6.8), linetype="longdash", size=1)+
  labs(title="", 
       x="PC2",
       fill="stage") +
  scale_fill_manual(values = c("Path1" = alpha("#662506",0.2),
                               'Path2'=alpha("#006d2c",0.2),
                               'Path3'= alpha("#08519c",0.2)),
                    breaks = c("Path1","Path2",'Path3'))+
  mytheme
dev.off()

pdf('./PATH.PC3.pdf')
g <- ggplot(as.data.frame(sim2_PCA_cor), aes(PC3))
g + geom_density(aes(fill=factor(layout$stages)), alpha=0.5) + 
  #geom_vline(xintercept=c(-6.8,-6.8), linetype="longdash", size=1)+
  labs(title="", 
       x="PC3",
       fill="stage") +
  scale_fill_manual(values = c("Path1" = alpha("#662506",0.2),
                               'Path2'=alpha("#006d2c",0.2),
                               'Path3'= alpha("#08519c",0.2)),
                    breaks = c("Path1","Path2",'Path3'))+
  mytheme
dev.off()

pdf('./PATH.PC4.pdf')
g <- ggplot(as.data.frame(sim2_PCA_cor), aes(PC4))
g + geom_density(aes(fill=factor(layout$stages)), alpha=0.5) + 
  #geom_vline(xintercept=c(-6.8,-6.8), linetype="longdash", size=1)+
  labs(title="", 
       x="PC4",
       fill="stage") +
  scale_fill_manual(values = c("Path1" = '#984ea3',
                               'Path2'='#377eb8',
                               'Path3'='#4daf4a'),
                    breaks = c("Path1","Path2",'Path3'))+
  mytheme
dev.off()

pdf('./PATH.PC5.pdf')
g <- ggplot(as.data.frame(sim2_PCA_cor), aes(PC5))
g + geom_density(aes(fill=factor(layout$stages)), alpha=0.5) + 
  #geom_vline(xintercept=c(-6.8,-6.8), linetype="longdash", size=1)+
  labs(title="", 
       x="PC5",
       fill="stage") +
  scale_fill_manual(values = c("Path1" = '#984ea3',
                               'Path2'='#377eb8',
                               'Path3'='#4daf4a'),
                    breaks = c("Path1","Path2",'Path3'))+
  mytheme
dev.off()
```


```R
labels <-  fetchLabels(MarkovObject=MarkovHC_PATH,
                       MarkovLevels=1:length(MarkovHC_PATH$hierarchicalStructure))

basins <- labels[,17]

for(i in 1:length(basins)){
    basins[i] <- str_split( basins[i], '\\+')[[1]][1]
}

basins <- paste('basin',basin,sep = '')
level <- 17
for (i in 1:length(MarkovHC_PATH$hierarchicalStructure[[level]]$attractorPoints)) {
  basins[MarkovHC_PATH$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
```


```R
layout <- sim2_plot$data%>%as.data.frame()

label <- findTransitionPoints(MarkovObject = MarkovHC_PATH,
                              level = level,
                              basinA = 1,
                              basinB = 2)
basins[which(label==1)] <- 'transition points\nbetween 1 and 2'

label <- findTransitionPoints(MarkovObject = MarkovHC_PATH,
                              level = level,
                              basinA = 3,
                              basinB = 1)
basins[which(label==1)] <- 'transition points\nbetween 3 and 1'

layout$basins <- basins
pdf(file = './PATH.B.pdf', width = 9, height = 5.5)
ggplot(data=layout, mapping =  aes(x=X, y=Y)) +
  geom_point(size=3, shape=21, aes(fill=basins), color="#525252")+
  xlim(min(layout$X)-1,max(layout$X)+1)+
  ylim(min(layout$Y)-1,max(layout$Y)+1)+
  mytheme+
  scale_fill_manual(
    values = c("attractors"=alpha("#e41a1c",1),
               "level17_basin1"=alpha("#377eb8",0.7),
               "level17_basin2"=alpha("#4daf4a",0.7),
               "level17_basin3"=alpha("#984ea3",0.7),
               "transition points\nbetween 1 and 2"=alpha("#ff7f00",1),
               "transition points\nbetween 3 and 1"=alpha("#ff7f00",1)),
    breaks = c("attractors",
               "level17_basin1",
               "level17_basin2",
               "level17_basin3",
               'transition points\nbetween 1 and 2',
               "transition points\nbetween 3 and 1"))
dev.off()
```


```R
save.image('./TransitionAndPath.RData')
```


```R
centrality_scores <- MarkovHC_PATH$midResults$centrality_scores
pdf(file = './PATH.centrality_scores.pdf', width = 7.5, height = 5.5)
ggplot(data=layout, mapping =  aes(x=X, y=Y, color=centrality_scores)) +
  geom_point(size=1, shape=19)+
  xlim(min(layout$X)-1,max(layout$X)+1)+
  ylim(min(layout$Y)-1,max(layout$Y)+1)+
  mytheme+
  xlab("X") + ylab("Y")
dev.off()
```

# Transition points and the path


```R
basin <- vector(length = nrow(MarkovHC_PATH$midResults$symmetric_KNN_graph))
level <- 17
MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_PATH,
                                   level = level,
                                   basinA = 3,
                                   basinB = 2)
#path points
Pathpoint <- c()
for(i in MarkovHCPath[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_PATH$hierarchicalStructure[[1]]$basinPoints[[i]])
}

label <- MarkovHCPath[[1]]
basin[which(label==1)] <- 'transition path\nfrom 3 to 2'

#find the transition points
for(i in 1:length(MarkovHCPath[[2]])){
  if((MarkovHCPath[[2]][i] %in% MarkovHC_PATH$hierarchicalStructure[[level]]$graphvertex_basins[[3]])&(!(MarkovHCPath[[2]][i+1] %in% MarkovHC_PATH$hierarchicalStructure[[level]]$graphvertex_basins[[3]]))){
    transitionPoint <- c(MarkovHC_PATH$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath[[2]][i]]],
                         MarkovHC_PATH$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath[[2]][i+1]]])
  }
}

PATHscaledData <- GetAssayData(object = PATHobject, slot = "scale.data")%>%t()
PATHscaledData <- subset(PATHscaledData, rownames(PATHscaledData)%in%names(Pathpoint))
PATHscaledData <- PATHscaledData[order(factor(rownames(PATHscaledData), levels = names(Pathpoint)), decreasing = FALSE),]
PATHscaledData <- as.data.frame(PATHscaledData)
plot(PATHscaledData[,which(colnames(PATHscaledData)=='Gene1')])
PATHscaledData$pseudotime <- 1:nrow(PATHscaledData)

transitionPoint_plot <- subset(PATHscaledData, rownames(PATHscaledData)%in%names(transitionPoint))

#get counts
PATHcounts <- GetAssayData(object = PATHobject, slot = "counts")%>%as.matrix()%>%t()

PATHcounts <- subset(PATHcounts, rownames(PATHcounts)%in%names(Pathpoint))
PATHcounts <- PATHcounts[order(factor(rownames(PATHcounts), levels = names(Pathpoint)), decreasing = FALSE),]
PATHcounts <- as.data.frame(PATHcounts)
plot(PATHcounts[,which(colnames(PATHcounts)=='Gene2')])
pseudotime <- as.data.frame(1:nrow(PATHcounts))
```

# use monocle to find differentially expressed genes along the path


```R
gene_metadata <- as.data.frame(colnames(PATHcounts))
rownames(gene_metadata) <- gene_metadata[,1]
gene_metadata $ gene_short_name <- gene_metadata[,1] 
colnames(gene_metadata) <- c('gene_short_name','ensembleID')
rownames(pseudotime) <- rownames(PATHcounts)

path_object <- newCellDataSet(as.matrix(t(PATHcounts)),
                              phenoData =  new("AnnotatedDataFrame", data = pseudotime),
                              featureData =  new("AnnotatedDataFrame", data = gene_metadata),
                              lowerDetectionLimit = 0.1,
                              expressionFamily = negbinomial.size())

path_object <- estimateSizeFactors(path_object)
path_object <- estimateDispersions(path_object)
path_object <- detectGenes(path_object, min_expr = 0)
path_pdata <- pData(path_object)
path_fdata <- fData(path_object)
path_expressed_genes <- row.names(subset(fData(path_object),num_cells_expressed >= 0))

colnames(pData(path_object))[1] <- 'Pseudotime'  
pData(path_object)$basins <- mapvalues(x = rownames(pData(path_object)), 
                                       from = rownames(PATHobject@meta.data), 
                                       to = PATHobject@meta.data$Group, 
                                       warn_missing = FALSE) 
path_object_subset <- path_object[c('Gene1'),]
plot_genes_in_pseudotime(path_object_subset, color_by = "basins",
                         cell_size = 2,relative_expr=TRUE)+ scale_color_manual(
  breaks = c("Path1", "Path2", "Path3"), 
  values=c('Path1'='#e41a1c', "Path2"='#377eb8', "Path3"='#4daf4a')) + theme(legend.position = "right")+
  mytheme

path_diff_test_res <- differentialGeneTest(path_object[PATHobject@assays$RNA@var.features,],
                                           fullModelFormulaStr = "~sm.ns(Pseudotime)")
path_sig_gene_names <- row.names(subset(path_diff_test_res, qval < 0.1))

path_pheatmap <- plot_pseudotime_heatmap(path_object[path_sig_gene_names,],
                                          num_clusters = 3,
                                          cores = 1,
                                          show_rownames = F,
                                          return_heatmap = TRUE)
plot_pseudotime_heatmap(path_object[path_pheatmap$tree_row$labels[path_pheatmap$tree_row[['order']]],],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = F,
                        cluster_rows=FALSE)

geneclusters <- cutree(path_pheatmap$tree_row, k=3)
orderedGenes <- c()
for(i in c(1,2,3)){
  orderedGenes <- c(orderedGenes, names(geneclusters[which(geneclusters==i)]))
}
pdf("./SplatterpseudotimeDiffGenes.pdf")
plot_pseudotime_heatmap(path_object[orderedGenes,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = F,
                        cluster_rows=FALSE)
dev.off()
save.image('./SplatterpseudotimeDiffGenes.RData')
#Gene
#Gene255
#Gene1901
#Gene2059

pdf('./Gene255.pdf', width = 3.5, height = 2)
ggplot(PATHscaledData,
       aes(x=pseudotime,y=Gene255)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=Gene255), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("Gene255")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

pdf('./Gene1901.pdf', width = 3.5, height = 2)
ggplot(PATHscaledData,
       aes(x=pseudotime,y=Gene1901)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=Gene1901), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("Gene1901")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

pdf('./Gene2059.pdf', width = 3.5, height = 2)
ggplot(PATHscaledData,
       aes(x=pseudotime,y=Gene2059)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=Gene2059), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("Gene2059")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()
```

# 3D visualization


```R
layout <- sim2_plot$data%>%as.data.frame()
layout$stages <- layout$colour_by
PCA123 <- Embeddings(object = PATHobject, reduction = "pca")[,1:3]
pdf('./splatter.groundtruth.pdf')
scatter3D(x=PCA123[,1],
          y=PCA123[,2], 
          z=PCA123[,3], 
          colvar = as.numeric(as.factor(layout$stages)), 
          col = c(alpha("#662506",0.2),alpha("#006d2c",0.2),  alpha("#08519c",0.2)),
          pch = 19, cex = 1.5,theta = 40, phi = 50)
dev.off()

#MarkovHC basins
labels <-  fetchLabels(MarkovObject=MarkovHC_PATH,
                       MarkovLevels=1:length(MarkovHC_PATH$hierarchicalStructure))

basins <- labels[,17]

for(i in 1:length(basins)){
    basins[i] <- str_split( basins[i], '\\+')[[1]][1]
}

basins <- paste('basin',basin,sep = '')
level <- 17
for (i in 1:length(MarkovHC_PATH$hierarchicalStructure[[level]]$attractorPoints)) {
  basins[MarkovHC_PATH$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
```


```R
pdf('./splatter.basins.pdf')
scatter3D(x=PCA123[,1],
          y=PCA123[,2], 
          z=PCA123[,3], 
          colvar = as.numeric(as.factor(basins)), 
          col = c(alpha("#ce1256",0.8),alpha("#006d2c",0.2),  alpha("#08519c",0.2),alpha("#662506",0.2)),
          pch = 19, cex = 1.5,theta = 40, phi = 50)
dev.off()

#PATH
basins <- character(length = ncol(mat2))
for (i in 1:length(basins)) {
  indexi <- subset(sankeyResult_PATH,sankeyResult_PATH[,1]==paste('index_basin',i,sep = ""))
  basins[i] <- paste(unique(indexi[,3])%>%sort(), collapse = '+')
}

label <- findTransitionPath(MarkovObject = MarkovHC_PATH,
                            level = level,
                            basinA = 3,
                            basinB = 2)
basins[which(label[[1]]==1)] <- 'transition path'


#find the transition points
transitionPoint <- c()
for(i in 1:length(label[[2]])){
  if((label[[2]][i] %in% MarkovHC_PATH$hierarchicalStructure[[level]]$graphvertex_basins[[3]])&(!(label[[2]][i+1] %in% MarkovHC_PATH$hierarchicalStructure[[level]]$graphvertex_basins[[3]]))){
    transitionPoint <- c(MarkovHC_PATH$hierarchicalStructure[[1]]$basinPoints[[label[[2]][i]]],
                         MarkovHC_PATH$hierarchicalStructure[[1]]$basinPoints[[label[[2]][i+1]]])
  }
}

basins[transitionPoint] <- 'critical_point'

pdf('./splatter.path3to2.pdf')
#x从小到大，basin3,basin1,basin2
scatter3D(x=PCA123[,1],
          y=PCA123[,2], 
          z=PCA123[,3], 
          colvar = as.numeric(as.factor(basins)), 
          col = c(#alpha("#ce1256",0.8),
                  alpha('#54278f',1),
                  alpha("#737373",0.1),
                  alpha("#737373",0.1),
                  alpha("#737373",0.1),
                  alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 50)
dev.off()

#density
#centrality scores
centrality_scores <- MarkovHC_PATH$midResults$centrality_scores
centrality_colors <- colorRampPalette(colors = c("#d9d9d9","#969696","#000000"))(length(centrality_scores))

pdf('./splatter.centrality.pdf')
scatter3D(x=PCA123[,1],
          y=PCA123[,2], 
          z=PCA123[,3], 
          colvar = centrality_scores, 
          col = centrality_colors,
          pch = 19, cex = 1.5,theta = 40, phi = 50)
dev.off()
```
