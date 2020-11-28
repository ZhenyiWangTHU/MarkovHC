```R
#Petropoulos-------------------------------------------------------------------
library(MarkovHC)
library(ggplot2)
library(Seurat)
library(ggalluvial)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MASS)
library(stringr)
library(pheatmap)
library(phateR)
library(plot3D)

petro <- read.table(file = "./Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos/E-MTAB-3929.processed.1/counts.txt", header = T)

mytheme <-  theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 1,
                                           colour = "black"),
                  axis.title.x =element_text(size = 20,
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
                  legend.title = element_blank())
#arrange the label-------------------------------------------------------------
petro_rowname  <- colnames(petro)
petro_rowname_label <- data.frame(1:1529, stringsAsFactors = F)
for (i in 1:length(petro_rowname)) {
  if(grepl("E3\\.[0-9]", petro_rowname[i])){petro_rowname_label[i,1] <- "E3"
  }else if(grepl("E4\\.[0-9]", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E4"
  }else if(grepl("E4\\.late", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E4 late"
  }else if(grepl("E5\\.early", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E5 early"
  }else if(grepl("E5\\.[0-9]", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E5"
  }else if(grepl("E6\\.[0-9]", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E6"
  }else if(grepl("E7\\.[0-9]", petro_rowname[i])){
    petro_rowname_label[i,1] <- "E7"
  }
}

#Seurat
ESCobject <- CreateSeuratObject(counts = petro,
                                 project = 'ESCobject',
                                 min.cells = 10,
                                 min.feature = 200)
ESCobject[["percent.mt"]] <- PercentageFeatureSet(ESCobject, pattern = "^mt-")
VlnPlot(ESCobject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ESCobjectplot1 <- FeatureScatter(ESCobject, feature1 = "nCount_RNA", feature2 = "percent.mt")
ESCobjectplot2 <- FeatureScatter(ESCobject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(ESCobjectplot1, ESCobjectplot2))

ESCobject <- NormalizeData(ESCobject, normalization.method = "LogNormalize", scale.factor = 10000)

ESCobject <- FindVariableFeatures(ESCobject, selection.method = "vst", nfeatures = 3000)
# Identify the 10 most highly variable genes
ESCobjecttop10 <- head(VariableFeatures(ESCobject), 10)
# plot variable features with and without labels
ESCobjectplot1 <- VariableFeaturePlot(ESCobject)
ESCobjectplot2 <- LabelPoints(plot = ESCobjectplot1, points = ESCobjecttop10, repel = TRUE)
ESCobjectplot2
ESCobject <- ScaleData(ESCobject, features = rownames(ESCobject))
ESCobject <- RunPCA(ESCobject, features = VariableFeatures(object = ESCobject), verbose=FALSE)

ElbowPlot(ESCobject, ndims = 50)

ESCobject <- RunTSNE(object = ESCobject, dims=1:20)
ESCobject <- RunUMAP(object = ESCobject, dims=1:20, n.neighbors=100)
DimPlot(ESCobject, reduction = "tsne", group.by = 'stage')
DimPlot(ESCobject, reduction = "pca", group.by = 'stage')
DimPlot(ESCobject, reduction = "umap", group.by = 'stage')

phate.ESCobject <- phate(subset(GetAssayData(object = ESCobject, slot = "scale.data"),
                        rownames(GetAssayData(object = ESCobject, slot = "scale.data"))%in%ESCobject@assays$RNA@var.features)%>%t(),
                              knn = 50,
                              npca=20,
                              t=5,
                           ndim=3)

plot(phate.ESCobject$embedding)

ESCobject@meta.data$stage <- petro_rowname_label[,1]
colorvector <- c( "E3"="#66c2a5",
                  "E4"="#fc8d62",
                  "E4 late"="#8da0cb",
                  "E5 early"="#e78ac3",
                  "E5"="#a6d854",
                  "E6"="#ffd92f",
                  "E7"="#1f78b4")

DimPlot(ESCobject, reduction = "pca",group.by='stage', cols = colorvector, pt.size=3,label = T,
        label.size = 6)

mytheme <- theme(panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(size = 1,
                                          colour = "black"),
                 axis.title.x =element_text(size=20),
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
                                             face = "bold"))

layout <- as.data.frame(phate.ESCobject$embedding)
layout$stage <- ESCobject@meta.data[,5]
pdf(file = './groudtruth.pdf', width = 3.5, height = 3.5)
ggplot(data=layout, mapping =  aes(x=PHATE1, y=PHATE2)) +
  geom_point(size=1, shape=21, aes(fill=stage), color=alpha("#525252",0))+
  xlim(min(layout$PHATE1)-1,max(layout$PHATE1)+1)+
  ylim(min(layout$PHATE2)-1,max(layout$PHATE2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c( "E3"=alpha("#DC0000FF",0.7),
                "E4"=alpha("#08519c",0.7),
                "E4 late"=alpha("#9ecae1",0.7),
                "E5 early"=alpha("#7E6148FF",0.7),
                "E5"=alpha("#C6A87F",0.7),
                "E6"=alpha("#00A087FF",0.7),
                "E7"=alpha("#F39B7FFF",0.7)),
    breaks = c( "E3",
                "E4",
                "E4 late",
                "E5 early",
                "E5",
                "E6",
                "E7"))
dev.off()

scales::show_col(c(alpha("#9ecae1",0.7),alpha("#DC0000FF",0.7),
                   alpha("#4DBBD5FF",0.7),alpha("#00A087FF",0.7),
                   alpha("#08519c",0.7),alpha("#7E6148FF",0.7),
                   alpha("#C6A87F",0.7),alpha("#F39B7FFF",0.7)))

#3D
pdf('./groundtruth.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = as.numeric(as.factor(layout$stage)), 
          col = c( alpha("#DC0000FF",0.7),
                   alpha("#08519c",0.7),
                   alpha("#9ecae1",0.7),
                   alpha("#7E6148FF",0.7),
                   alpha("#C6A87F",0.7),
                   alpha("#00A087FF",0.7),
                   alpha("#F39B7FFF",0.7)),
          pch = 19, 
          cex = 0.5,
          theta = 0,
          phi = 125)
dev.off()
#MarkovHC------------------------------------------------------------------
scale.data <- GetAssayData(object = ESCobject, slot = "scale.data")
scale.data <- subset(scale.data, rownames(scale.data)%in%ESCobject@assays$RNA@var.features)

MarkovHC_ESCobject <- MarkovHC(origin_matrix=Embeddings(object = ESCobject, reduction = "pca")[,1:20]%>%t(),
                                transformtype="none",
                                KNN=50,
                                basecluster="kmeans",
                                dobasecluster=TRUE,
                                baseclusternum=500,
                                emphasizedistance=1,
                                weightDist=2,
                                weightDens=0.5,
                                cutpoint=0.01,
                                showprocess=FALSE,
                                bn=2,
                                minBasinSize=0.2,
                                noiseBasinSize=10)
save(MarkovHC_ESCobject,file = './MarkovHC_ESCobject.Rdata')

#Figure3.B
labels <-  fetchLabels(MarkovObject=MarkovHC_ESCobject,
                       MarkovLevels=1:15)

#level6
layout <- as.data.frame(phate.ESCobject$embedding)
layout$basins <- labels[,6]
for(i in 1:nrow(layout)){
  layout$basins[i] <- str_split(layout$basins[i],'\\+')[[1]][1]
}
pdf(file = './level6withlegend.pdf', width = 3.5, height = 3.5)
 ggplot(data=layout, mapping =  aes(x=PHATE1, y=PHATE2)) +
  geom_point(size=2, shape=21, aes(fill=basins), color=alpha("#525252",0))+
  xlim(min(layout$PHATE1)-1,max(layout$PHATE1)+1)+
  ylim(min(layout$PHATE2)-1,max(layout$PHATE2)+1)+
  mytheme+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values =c("3" = alpha("#DC0000FF",0.7),
              "1" = alpha("#08519c",0.7),
              "7" = alpha("#9ecae1",0.7),
              "4" = alpha("#7E6148FF",0.7),
              "5" = alpha("#C6A87F",0.7),
              "9" = alpha("#00A087FF",0.7),
              "10" = alpha("#F39B7FFF",0.7),
              "2" = alpha('#6a51a3',0.7),
              "8" = alpha('#c994c7',0.7),
              "6" = alpha('#006d2c',0.7)),
    breaks = c("1",
               "2",
               "3",
               "4",
               "5",
               "6",
               "7",
               "8",
               "9",
               "10"))
dev.off()
#3D
colors <- c( alpha("#DC0000FF",0.7),
             alpha("#08519c",0.7),
             alpha("#4292c6",0.7),
             alpha("#7E6148FF",0.7),
             alpha("#bf812d",0.7),
             alpha("#00A087FF",0.7),
             alpha("#F39B7FFF",0.7),
             alpha('#6a51a3',0.7),
             alpha('#c994c7',0.7),
             alpha('#006d2c',0.7))



pdf('./level6.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = as.numeric(factor(labels[,6],levels=c(3,1,7,4,5,9,10,2,8,6))), 
          col = colors,
          pch = 19, 
          cex = 0.5,
          theta = 0,
          phi = 125)
dev.off()

#level4
for(j in 1:ncol(labels)){
  for(i in 1:nrow(labels)){
  labels[i,j] <- str_split(labels[i,j],'\\+')[[1]][1]
  }
}

dm_temp <- dist(Embeddings(object = ESCobject, reduction = "pca")[,1:20],method = "minkowski",p=2)
dm_matrix_temp <- as.matrix(dm_temp)

NeatLabels <- fetchNeatLabels(MarkovObject=MarkovHC_ESCobject,
                             level1=4,
                             level2=6,
                             labelDataFrame=labels,
                             dm_matrix=dm_matrix_temp,
                             minimumBasin=20)

NeatLabels[which(!(NeatLabels%in%c('9','12','16')))] <- 0

pdf('./ICMonlevel4.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = as.numeric(factor(NeatLabels)), 
          col = c(alpha("#737373",0.1),
                  alpha("#fb8072",1),
                  alpha('#c994c7',1),
                  alpha('#8dd3c7',1)),
          pch = 19, 
          cex = 0.5,
          theta = 0,
          phi = 125)
dev.off()
#tree
#plot the tree
pdf(file= './treeuseFunction.pdf', width = 15, height = 10)
plotHierarchicalStructure(MarkovObject=MarkovHC_ESCobject,
                          MarkovLevels=6:15,
                          colorVector=c( alpha("#08519c",0.7), 
                                         alpha("#6a51a3",0.7), 
                                         alpha("#DC0000FF",0.7),
                                         alpha("#7E6148FF",0.7),
                                         alpha("#bf812d",0.7), 
                                         alpha('#006d2c',0.7),
                                         alpha('#4292c6',0.7),  
                                         alpha("#c994c7",0.7),
                                         alpha('#00A087FF',0.7), 
                                         alpha("#F39B7FFF",0.7)
                          ))
dev.off()

#heatmap
#Find differentially expressed genes
layout$basins <- labels[,6]
for(i in 1:nrow(layout)){
  layout$basins[i] <- str_split(layout$basins[i],'\\+')[[1]][1]
}
Idents(object = ESCobject) <- layout$basins

ESCobject.markers <- FindAllMarkers(ESCobject,
                                    min.pct = 0.25,
                                    logfc.threshold = 0.25,
                                    only.pos=TRUE)
ESCobject.markerstop50 <- ESCobject.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>%as.data.frame()

ordergenes_expression_matrix <- GetAssayData(object = ESCobject, slot = "scale.data")%>%as.data.frame()
ESCobject.markerstop50 <- ESCobject.markerstop50[order(factor(ESCobject.markerstop50$cluster, levels = c("3","1","7","5","4","8","2","6","9","10")), decreasing = FALSE),]
ESCobject.markerstop50$cluster <- plyr::mapvalues(ESCobject.markerstop50$cluster, 
                                                       from=c("3","1","7","5","4","8","2","6","9","10"),
                                                       to=c('basin1',
                                                            'basin2',
                                                            'basin3',
                                                            'basin4',
                                                            'basin5',
                                                            'basin6',
                                                            'basin7',
                                                            'basin8',
                                                            'basin9',
                                                            'basin10'))
write.table(ESCobject.markerstop50, file = './Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos/markerGenes.txt', quote = FALSE, sep=' ', row.names = F)


ordergenes_expression_matrix <- subset(ordergenes_expression_matrix, rownames(ordergenes_expression_matrix)%in%ESCobject.markerstop50$gene)
ordergenes_expression_matrix <- ordergenes_expression_matrix[order(factor(rownames(ordergenes_expression_matrix), levels = unique(ESCobject.markerstop50$gene))),]

annotation_col_C_andcluster = data.frame(Days=factor(ESCobject@meta.data$stage),
                                         Basins=factor(Idents(ESCobject)))
rownames(annotation_col_C_andcluster) = colnames(ordergenes_expression_matrix)
ann_colors_C = list(
  Days =  c("E3"="#DC0000FF",
  "E4"="#08519c",
  "E4 late"="#9ecae1",
  "E5 early"="#7E6148FF",
  "E5"="#C6A87F",
  "E6"="#00A087FF",
  "E7"="#F39B7FFF"),
  Basins = c("1"="#08519c",
             "2"="#6a51a3",
             "3"="#DC0000FF",
             "4"="#7E6148FF",
             "5"="#bf812d",
             "6"="#006d2c",
             "7"="#4292c6",
             "8"="#c994c7",
             "9"="#00A087FF",
             "10"="#F39B7FFF"))

ordered_genes_expression_matrix <- t(ordergenes_expression_matrix)%>%as.data.frame()
ordered_genes_expression_matrix$Basins <- factor(Idents(ESCobject), levels = rev(c("3","1","7","5","4","8","2","6","9","10")))
ordered_genes_expression_matrix$Days <- factor(ESCobject@meta.data$stage, levels = rev(c( "E3",
                                                                                      "E4",
                                                                                      "E4 late",
                                                                                      "E5 early",
                                                                                      "E5",
                                                                                      "E6",
                                                                                      "E7")))
ordered_genes_expression_matrix <- doBy::orderBy(~ Basins + Days, ordered_genes_expression_matrix)
ordered_genes_expression_matrix <- ordered_genes_expression_matrix[,-which(colnames(ordered_genes_expression_matrix)%in%c('Basins','Days'))]
ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)
ordered_genes_expression_matrix_copy <- ordered_genes_expression_matrix

ordered_genes_expression_matrix[ordered_genes_expression_matrix>4] <- 4
ordered_genes_expression_matrix[ordered_genes_expression_matrix< (-4)] <- (-4)
ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)
ordered_genes_expression_matrix <- scale(ordered_genes_expression_matrix, center=T)
ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)
pheatmap(as.matrix(ordered_genes_expression_matrix), cluster_rows = F, cluster_cols =F,
         scale = "none" ,
         legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                    max(ordered_genes_expression_matrix),0.01)),
         color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
         breaks= seq(min(ordered_genes_expression_matrix),
                     max(ordered_genes_expression_matrix),
                     by=0.01),
         show_colnames = F, show_rownames = F,
         annotation_col  = annotation_col_C_andcluster,
         annotation_colors = ann_colors_C,width=15, heigheight = 10,
         fontsize =15,
         filename = './markerGeneHeatmap.pdf'
         )
# barplot
barplot_dataframe <- as.data.frame(colnames(ordered_genes_expression_matrix))
barplot_dataframe$basins <- mapvalues(colnames(ordered_genes_expression_matrix),
                                      from=rownames(annotation_col_C_andcluster),
                                      to= as.character(annotation_col_C_andcluster$Basins))
barplot_dataframe$days <- mapvalues(colnames(ordered_genes_expression_matrix),
                                      from=rownames(annotation_col_C_andcluster),
                                      to= as.character(annotation_col_C_andcluster$Days))
head(barplot_dataframe)
ggplot(barplot_dataframe,aes(x=days)) + 
  geom_bar(aes(fill=factor(basins)),position="fill")+scale_fill_manual(values=c("1"="#08519c",
                                                                                "2"="#6a51a3",
                                                                                "3"="#DC0000FF",
                                                                                "4"="#7E6148FF",
                                                                                "5"="#bf812d",
                                                                                "6"="#006d2c",
                                                                                "7"="#4292c6",
                                                                                "8"="#c994c7",
                                                                                "9"="#00A087FF",
                                                                                "10"="#F39B7FFF"))+
  mytheme+ylab('percentage')+guides(fill=guide_legend(title='phase'))+ theme(axis.text.x = element_text(size = 20,
                                                                                                        family = "sans",
                                                                                                        color = "black",
                                                                                                        face = "bold",
                                                                                                        vjust = 1,
                                                                                                        hjust = 1,
                                                                                                        angle = 45))


#GO
ESCobject.markerstop50$cluster <- as.character(ESCobject.markerstop50$cluster)
basins <- unique(ESCobject.markerstop50[,6])
for(i in 1:length(basins)){
  upregulatedGenes <- (subset(ESCobject.markerstop50, ESCobject.markerstop50[,6]==basins[i])%>%as.data.frame())[,7]
  GO_upregulatedGenes <- enrichGO(gene = upregulatedGenes,
                                  keyType = "SYMBOL",
                                  OrgDb = 'org.Hs.eg.db',
                                  ont = "BP",
                                  pAdjustMethod = "fdr",
                                  pvalueCutoff = 0.05,
                                  qvalueCutoff  = 0.2,
                                  minGSSize = 3,
                                  maxGSSize = 500,
                                  readable = FALSE)
  GO_upregulatedGenes.result <- as.data.frame(GO_upregulatedGenes@result)
  write.table(GO_upregulatedGenes.result, file = paste('./Single-Cell RNA-Seq Reveals Lineage and X Chromosome Dynamics in Human Preimplantation Embryos/',as.character(basins[i]),'.txt', sep=''))
}

#centrality scores
centrality_colors <- colorRampPalette(colors = c("#d9d9d9","#969696","#000000"))(1529)
centrality_scores <- MarkovHC_ESCobject$midResults$centrality_scores

pdf(file = './centrality_scores.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = centrality_scores, 
          col = centrality_colors,
          pch = 19, 
          cex = 0.5,
          theta = 0,
          phi = 125)

dev.off()

save.image('./ESCDevelopment/ESCDevelopment.RData')


#Seurat plot
ESCobject@meta.data$basin <- sankeyResult_ESCobject[,4]

DimPlot(ESCobject, reduction = "pca", group.by='basin', cols = colorvector, pt.size=3,label = T,
        label.size = 6)

pdf('./ESCDevelopment/DoHeatmap.pdf')
print(DoHeatmap(ESCobject, features = ESCobject.markerstop50$gene,angle = 30,size=3))
dev.off()

##transfer levels17 basin7 to that on level 15
basins_level17Andlevel15 <- sankeyResult_ESCobject[,4]
basins_level17Andlevel15[which(basins_level17Andlevel15=='level17_basin7')] <- subset(sankeyResult_ESCobject, sankeyResult_ESCobject[,4]=='level17_basin7')$level15
ESCobject@meta.data$bestBasins <- basins_level17Andlevel15
cellsOfLevel17Basin7 <- subset(x = ESCobject, subset = bestBasins%in%c("level15_basin7", "level15_basin8", "level15_basin12"))

Idents(object = cellsOfLevel17Basin7) <- cellsOfLevel17Basin7@meta.data$bestBasins

cellsOfLevel17Basin7.markers <- FindAllMarkers(cellsOfLevel17Basin7,
                                              #features=ESCobject@assays$RNA@var.features,
                                              min.pct = 0.25,
                                              logfc.threshold = 0.25,
                                              only.pos=TRUE)
cellsOfLevel17Basin7.markerstop50 <- cellsOfLevel17Basin7.markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)

cellsOfLevel17Basin7_matrix <- GetAssayData(object = cellsOfLevel17Basin7, slot = "scale.data")
cellsOfLevel17Basin7.markerstop50 <- cellsOfLevel17Basin7.markerstop50[order(factor(cellsOfLevel17Basin7.markerstop50$cluster, 
                                                                                    levels = c("level15_basin7",
                                                                                               "level15_basin8",
                                                                                               "level15_basin12"))),]

cellsOfLevel17Basin7_matrix <- subset(cellsOfLevel17Basin7_matrix, rownames(cellsOfLevel17Basin7_matrix)%in%cellsOfLevel17Basin7.markerstop50$gene)
cellsOfLevel17Basin7_matrix <- cellsOfLevel17Basin7_matrix[order(factor(rownames(cellsOfLevel17Basin7_matrix), levels = unique(cellsOfLevel17Basin7.markerstop50$gene))),]

celltypecahracter <- as.character(rep('level17_basin7', ncol(cellsOfLevel17Basin7_matrix)))

annotation_col_C_andcluster = data.frame(level17=factor(celltypecahracter),
                                         level15=factor(cellsOfLevel17Basin7@meta.data$bestBasins))
rownames(annotation_col_C_andcluster) = colnames(cellsOfLevel17Basin7_matrix)
ann_colors_C = list(
  level17 =  c(level17_basin7="#a6d854"),
  level15 = c(level15_basin7="#b2abd2",
             level15_basin8="#8073ac",
             level15_basin12="#542788"))

unique_combine_basin <- c("level15_basin7",
                          "level15_basin8",
                          "level15_basin12")
cellsOfLevel17Basin7_matrix <- t(cellsOfLevel17Basin7_matrix)
ordered_genes_expression_matrix <- subset(cellsOfLevel17Basin7_matrix, cellsOfLevel17Basin7@meta.data$bestBasins%in%unique_combine_basin[1])
for (i in 2:length(unique_combine_basin)) {
  subcluster <- subset(cellsOfLevel17Basin7_matrix, cellsOfLevel17Basin7@meta.data$bestBasins%in%unique_combine_basin[i])
  ordered_genes_expression_matrix <- rbind(ordered_genes_expression_matrix, subcluster)
}
ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)

ordered_genes_expression_matrix_copy <- ordered_genes_expression_matrix

#log_ordered_genes_expression_matrix <- log(ordered_genes_expression_matrix+1)
#log_ordered_genes_expression_matrix[log_ordered_genes_expression_matrix>4] <- 4
#ordered_genes_expression_matrix <- log_ordered_genes_expression_matrix
range(ordered_genes_expression_matrix)
pheatmap(as.matrix(ordered_genes_expression_matrix), cluster_rows = F, cluster_cols =F,
         scale = "row" ,
         legend_breaks= ceiling(seq(-4,
                                    4,0.01)),
         color = colorRampPalette(colors = c("turquoise1","black","gold"))(length(seq(-4,4,0.01))),
         breaks= seq(-4,
                     4,
                     by=0.01),
         show_colnames = F, show_rownames = F,
         annotation_col  = annotation_col_C_andcluster,
         annotation_colors = ann_colors_C,
         fontsize =15, filename = './ESdevelopmentlevel15.heatmap.pdf')

upregulatedGenes <- (subset(cellsOfLevel17Basin7.markers, cellsOfLevel17Basin7.markers[,6]=="level15_basin12")%>%as.data.frame())[,7]
GO_upregulatedGenes <- enrichGO(gene = upregulatedGenes,
                                keyType = "SYMBOL",
                                OrgDb = 'org.Hs.eg.db',
                                ont = "BP",
                                pAdjustMethod = "fdr",
                                pvalueCutoff = 0.05,
                                qvalueCutoff  = 0.2,
                                minGSSize = 3,
                                maxGSSize = 500,
                                readable = FALSE)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
#tree plot
#plot the tree
pdf(file= './ESdevelopment.treeuseFunction.pdf', width = 15, height = 10)
plotHierarchicalStructure(MarkovObject=MarkovHC_ESCobject,
                          MarkovLevels=17:22,
                          colorVector=c("#fc8d62","#33a02c","#e31a1c","#ffd92f","#1f78b4",
                                        "#a6cee3","#a6d854","#66c2a5","#bf5b17","#e78ac3"))
dev.off()

#marker genes
ESCobject@meta.data$basin <- labels$level6
for(i in 1:nrow(ESCobject@meta.data)){
  ESCobject@meta.data$basin[i] <- str_split(ESCobject@meta.data$basin[i],'\\+')[[1]][1]
}
ESCobject@meta.data$basin <- factor(ESCobject@meta.data$basin, 
                                    levels = rev(c("3","1","7","5","8","4","2","6","9","10")))
ESCobject.markerstop3 <- ESCobject.markers %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC)%>%as.data.frame()
ESCobject.markerstop3 <- ESCobject.markerstop3[order(factor(ESCobject.markerstop3$cluster, levels = rev(c("3","1","7","5","4","8","2","6","9","10"))), decreasing = FALSE),]

pdf('./MarkerGenes.pdf', width = 5, height = 5)
DotPlot(ESCobject, group.by='basin', 
        features = unique(rev(c('SOX17','SOX7','GATA6',
                                'GATA2','GATA3','SOX2',
                                'PDGFRA',
                                #mural marekr genes
                                'SLC15A2', 'PDCD4','PLCE1','WIPF3',
                                #polar marker genes
                                'SLCO4A1', 'NRP1','GCM1')))) + coord_flip()+NoLegend()
dev.off()

#marker genes of ICM
ESCobject@meta.data$basin <- NeatLabels

ESCobject@meta.data$basin <- factor(ESCobject@meta.data$basin, 
                                    levels = rev(c("0","16","12","9")))
pdf('./MarkerGenesICM.pdf', width = 5, height = 5)
DotPlot(ESCobject, group.by='basin', 
        features = unique(rev(c('SOX2','TDGF1','DPPA5','GDF3','PRDM14',
                                "PDGFRA",'FGFR2','LAMA4','HNF1B')))) + coord_flip()#+NoLegend()
dev.off()

#The sequence of C_cut
pdf(file = './C_cut_seq1_11.pdf', width = 3.5, height = 3.5)
ggplot(C_cut_seq, aes(x=level, y=C_cut_seq))+ geom_point(size=1,shape=19)+mytheme
dev.off()

```

# path


```R
#find the transition paths
basin <- vector(length = nrow(MarkovHC_ESCobject$midResults$symmetric_KNN_graph))
level <- 6
MarkovHCPath31 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 3,
                                   basinB = 1)
MarkovHCPath17 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 1,
                                   basinB = 7)
MarkovHCPath75 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                     level = level,
                                     basinA = 7,
                                     basinB = 5)
MarkovHCPath58 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 5,
                                   basinB = 8)
#path length matrix
#level 6
pathLength <- matrix(0,10,10)
for(i in 1:10){
  for(j in 1:10){
    MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                       level = level,
                                       basinA = i,
                                       basinB = j)
    pathLength[i,j] <- MarkovHCPath[[3]]
  }
}
pathLength <- pathLength[c(3,1,7,5,4,8,2,6,9,10),
                         c(3,1,7,5,4,8,2,6,9,10)]
pheatmap::pheatmap(as.matrix(pathLength), cluster_rows = F, cluster_cols =F,
         scale = "none" ,
         legend_breaks= ceiling(seq(min(pathLength),
                                    max(pathLength),0.01)),
         color = colorRampPalette(colors = c("#e41a1c","#377eb8","#deebf7"))(length(seq(min(pathLength),max(pathLength),0.01))),
         breaks= seq(min(pathLength),
                     max(pathLength),
                     by=0.01),
         show_colnames = F, show_rownames = F,
         fontsize =10,
         cellwidth = 30,
         cellheight = 30,
         display_numbers=TRUE,
         number_color = 'black',
         number_format = "%.3f"
)

#level 3
pathLength <- matrix(0,42,42)
level <- 4
for(i in 1:42){
  for(j in 1:42){
    MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                       level = level,
                                       basinA = i,
                                       basinB = j)
    if(length(MarkovHCPath)>0){pathLength[i,j] <- MarkovHCPath[[3]]}else{
      pathLength[i,j] <- Inf
    }
    
  }
}
pathLength <- pathLength[c(16,9,12),
                         c(16,9,12)]
pheatmap::pheatmap(as.matrix(pathLength), cluster_rows = F, cluster_cols =F,
                   scale = "none" ,
                   legend_breaks= ceiling(seq(min(pathLength),
                                              1,0.01)),
                   color = colorRampPalette(colors = c("#e41a1c","#377eb8","#deebf7"))(length(seq(min(pathLength),1,0.01))),
                   breaks= seq(min(pathLength),
                               1,
                               by=0.01),
                   show_colnames = F, show_rownames = F,
                   width=10, 
                   heigheight = 10,
                   fontsize =10,
                   cellwidth = 30,
                   cellheight = 30,
                   display_numbers=TRUE,
                   number_color = 'black',
                   number_format = "%.3f",
                   filename = './HeatmapLvel4.pdf'
)

#path points
Pathpoint <- c()
for(i in MarkovHCPath31[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath17[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath75[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath58[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}

label31 <- MarkovHCPath31[[1]]
label17 <- MarkovHCPath17[[1]]
label75 <- MarkovHCPath75[[1]]
label58 <- MarkovHCPath58[[1]]
basin[((label31==1)|(label17==1)|(label75==1)|(label58==1))] <- 'transition path\nfrom 3 to 8'


#find the critical points
for(i in 1:length(MarkovHCPath31[[2]])){
  if((MarkovHCPath31[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[3]])&(!(MarkovHCPath31[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[3]]))){
    transitionPoint31 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath31[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath31[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath17[[2]])){
  if((MarkovHCPath17[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[1]])&(!(MarkovHCPath17[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[1]]))){
    transitionPoint17 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath17[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath17[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath75[[2]])){
  if((MarkovHCPath75[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[7]])&(!(MarkovHCPath75[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[7]]))){
    transitionPoint75 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath75[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath75[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath58[[2]])){
  if((MarkovHCPath58[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[5]])&(!(MarkovHCPath58[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[5]]))){
    transitionPoint58 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath58[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath58[[2]][i+1]]])
  }
}


basin[transitionPoint58] <- 'transitionPoint'
transitionPoint <- transitionPoint58

ESCscaledData <- GetAssayData(object = ESCobject, slot = "scale.data")%>%t()
ESCscaledData <- subset(ESCscaledData, rownames(ESCscaledData)%in%names(Pathpoint))
ESCscaledData <- ESCscaledData[order(factor(rownames(ESCscaledData), levels = unique(names(Pathpoint))), decreasing = FALSE),]
ESCscaledData <- as.data.frame(ESCscaledData)
plot(ESCscaledData[,which(colnames(ESCscaledData)=='PTGES')])
ESCscaledData$pseudotime <- 1:nrow(ESCscaledData)

transitionPoint_plot <- subset(ESCscaledData, rownames(ESCscaledData)%in%unique(names(transitionPoint)))
#ESCscaledData <- ESCscaledData[-which(rownames(ESCscaledData)%in%names(transitionPoint)),]

#SOX2
pdf('./ESpseudotimeSOX2.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=SOX2)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=SOX2), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("SOX2")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#NANOG
pdf('./ESpseudotimeNANOG.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=NANOG)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=NANOG), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("NANOG")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#GATA2
pdf('./ESpseudotimeGATA2.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=GATA2)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=GATA2), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("GATA2")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#TDGF1
TDGF1 <- ggplot(ESCscaledData,
       aes(x=pseudotime,y=TDGF1)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=TDGF1), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("TDGF1")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme



#selected Genes
selectedGenes <- c('FRMD4B','GATA3','PDGFRA','GATA6','DPPA5','FGF4',
                   'MYC','PTGES','LRP2','HNF1B','HNF4A','ARGFX','TDGF1',
                   'SLC34A2','PDGFA','WNT7A','APOA1','P4HA1','SOX2','NODAL',
                   'NANOG','GATA2')
for(i in selectedGenes){
  pdf(paste('./ESpseudotime',i,'.pdf', sep=''),width = 3.5, height = 2)
  print(ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=i)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=i), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme)
  dev.off()
}

selectedGenes <- c("GATA2","PDGFRA", 
                   "DPPA5", "MYC",
                   "PTGES","LRP2","HNF1B",
                   "ARGFX","TDGF1","SLC34A2",
                   "PDGFA","P4HA1")
for(i in selectedGenes){
  pdf(paste('./ESpseudotime',i,'.pdf', sep=''),width = 3.5, height = 2)
  print(ggplot(ESCscaledData,
               aes_string(x="pseudotime",y=i)) + geom_point(shape=19, size=1.5, color='#feb24c') +
          geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=i), color='#54278f')+
          ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
          geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme)
  dev.off()
}


selectedGenes <- c("GATA2","PDGFA","PTGES",
                   "MYC","SLC34A2","ARGFX")
x <- c(1,2,3,1,2,3)
y <- c(1,1,1,2,2,2)
pdf(file = './combine_allscatters.pdf', width = 9, height = 6)
grid.newpage()  ###新建图表版面
pushViewport(viewport(layout = grid.layout(3,2))) ####将版面分面
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
for(i in 1:length(selectedGenes)){
print(ggplot(ESCscaledData,
                   aes_string(x="pseudotime",y=selectedGenes[i])) + geom_point(shape=19, size=1.5, color='#feb24c') +
              geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[i]), color='#54278f')+
              ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[i])+
              geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme , vp = vplayout(x[i],y[i]))
}
dev.off()

p1 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[1])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[1]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[1])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p2 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[2])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[2]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[2])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p3 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[3])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[3]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[3])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p4 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[4])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[4]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[4])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p5 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[5])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[5]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[5])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p6 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[6])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[6]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[6])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
pdf(file = './combine_allscatters.pdf', width = 6, height = 9)
ggpubr::ggarrange(
  p1,p5,p3,p4,p2,p6,
  ncol = 2, nrow = 3,align = 'v')
dev.off()  
  
#plot the path
basin[which(basin==FALSE)] <- 0
basin[which(basin=='transition path\nfrom 3 to 8')] <- 1
basin[which(basin=='transitionPoint')] <- 2
basin <- as.numeric(basin)
pdf(file = './ESPCpathplot.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = basin, 
          col = c(alpha("#737373",0.1),
            alpha("#feb24c",1),
            alpha('#54278f',1)),
          pch = 19, 
          cex = 1,
          theta = 0,
          phi = 125)
dev.off()

save.image('./ESCDevelopmentPath.RData')

```

# path


```R
#find the transition paths
basin <- vector(length = nrow(MarkovHC_ESCobject$midResults$symmetric_KNN_graph))
level <- 6
MarkovHCPath31 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 3,
                                   basinB = 1)
MarkovHCPath17 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 1,
                                   basinB = 7)
MarkovHCPath75 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                     level = level,
                                     basinA = 7,
                                     basinB = 5)
MarkovHCPath58 <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                   level = level,
                                   basinA = 5,
                                   basinB = 8)
#path length matrix
#level 6
pathLength <- matrix(0,10,10)
for(i in 1:10){
  for(j in 1:10){
    MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                       level = level,
                                       basinA = i,
                                       basinB = j)
    pathLength[i,j] <- MarkovHCPath[[3]]
  }
}
pathLength <- pathLength[c(3,1,7,5,4,8,2,6,9,10),
                         c(3,1,7,5,4,8,2,6,9,10)]
pheatmap::pheatmap(as.matrix(pathLength), cluster_rows = F, cluster_cols =F,
         scale = "none" ,
         legend_breaks= ceiling(seq(min(pathLength),
                                    max(pathLength),0.01)),
         color = colorRampPalette(colors = c("#e41a1c","#377eb8","#deebf7"))(length(seq(min(pathLength),max(pathLength),0.01))),
         breaks= seq(min(pathLength),
                     max(pathLength),
                     by=0.01),
         show_colnames = F, show_rownames = F,
         fontsize =10,
         cellwidth = 30,
         cellheight = 30,
         display_numbers=TRUE,
         number_color = 'black',
         number_format = "%.3f"
)

#level 3
pathLength <- matrix(0,42,42)
level <- 4
for(i in 1:42){
  for(j in 1:42){
    MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_ESCobject,
                                       level = level,
                                       basinA = i,
                                       basinB = j)
    if(length(MarkovHCPath)>0){pathLength[i,j] <- MarkovHCPath[[3]]}else{
      pathLength[i,j] <- Inf
    }
    
  }
}
pathLength <- pathLength[c(16,9,12),
                         c(16,9,12)]
pheatmap::pheatmap(as.matrix(pathLength), cluster_rows = F, cluster_cols =F,
                   scale = "none" ,
                   legend_breaks= ceiling(seq(min(pathLength),
                                              1,0.01)),
                   color = colorRampPalette(colors = c("#e41a1c","#377eb8","#deebf7"))(length(seq(min(pathLength),1,0.01))),
                   breaks= seq(min(pathLength),
                               1,
                               by=0.01),
                   show_colnames = F, show_rownames = F,
                   width=10, 
                   heigheight = 10,
                   fontsize =10,
                   cellwidth = 30,
                   cellheight = 30,
                   display_numbers=TRUE,
                   number_color = 'black',
                   number_format = "%.3f"
)



#path points
Pathpoint <- c()
for(i in MarkovHCPath31[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath17[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath75[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}
for(i in MarkovHCPath58[[2]]){
  Pathpoint <- c(Pathpoint, MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[i]])
}

label31 <- MarkovHCPath31[[1]]
label17 <- MarkovHCPath17[[1]]
label75 <- MarkovHCPath75[[1]]
label58 <- MarkovHCPath58[[1]]
basin[((label31==1)|(label17==1)|(label75==1)|(label58==1))] <- 'transition path\nfrom 3 to 8'

#find the critical points
for(i in 1:length(MarkovHCPath31[[2]])){
  if((MarkovHCPath31[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[3]])&(!(MarkovHCPath31[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[3]]))){
    transitionPoint31 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath31[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath31[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath17[[2]])){
  if((MarkovHCPath17[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[1]])&(!(MarkovHCPath17[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[1]]))){
    transitionPoint17 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath17[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath17[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath75[[2]])){
  if((MarkovHCPath75[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[7]])&(!(MarkovHCPath75[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[7]]))){
    transitionPoint75 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath75[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath75[[2]][i+1]]])
  }
}

for(i in 1:length(MarkovHCPath58[[2]])){
  if((MarkovHCPath58[[2]][i] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[5]])&(!(MarkovHCPath58[[2]][i+1] %in% MarkovHC_ESCobject$hierarchicalStructure[[level]]$graphvertex_basins[[5]]))){
    transitionPoint58 <- c(MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath58[[2]][i]]],
                         MarkovHC_ESCobject$hierarchicalStructure[[1]]$basinPoints[[MarkovHCPath58[[2]][i+1]]])
  }
}

basin[transitionPoint58] <- 'transitionPoint'
transitionPoint <- transitionPoint58

ESCscaledData <- GetAssayData(object = ESCobject, slot = "scale.data")%>%t()
ESCscaledData <- subset(ESCscaledData, rownames(ESCscaledData)%in%names(Pathpoint))
ESCscaledData <- ESCscaledData[order(factor(rownames(ESCscaledData), levels = unique(names(Pathpoint))), decreasing = FALSE),]
ESCscaledData <- as.data.frame(ESCscaledData)
plot(ESCscaledData[,which(colnames(ESCscaledData)=='PTGES')])
ESCscaledData$pseudotime <- 1:nrow(ESCscaledData)

transitionPoint_plot <- subset(ESCscaledData, rownames(ESCscaledData)%in%unique(names(transitionPoint)))

#SOX2
pdf('./ESpseudotimeSOX2.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=SOX2)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=SOX2), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("SOX2")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#NANOG
pdf('./ESpseudotimeNANOG.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=NANOG)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=NANOG), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("NANOG")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#GATA2
pdf('./ESpseudotimeGATA2.pdf',width = 3.5, height = 2)
ggplot(ESCscaledData,
       aes(x=pseudotime,y=GATA2)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=GATA2), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("GATA2")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme
dev.off()

#TDGF1
TDGF1 <- ggplot(ESCscaledData,
       aes(x=pseudotime,y=TDGF1)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes(x=pseudotime,y=TDGF1), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle("TDGF1")+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme



#selected Genes
selectedGenes <- c('FRMD4B','GATA3','PDGFRA','GATA6','DPPA5','FGF4',
                   'MYC','PTGES','LRP2','HNF1B','HNF4A','ARGFX','TDGF1',
                   'SLC34A2','PDGFA','WNT7A','APOA1','P4HA1','SOX2','NODAL',
                   'NANOG','GATA2')
for(i in selectedGenes){
  pdf(paste('./ESpseudotime',i,'.pdf', sep=''),width = 3.5, height = 2)
  print(ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=i)) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=i), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme)
  dev.off()
}

selectedGenes <- c("GATA2","PDGFRA", 
                   "DPPA5", "MYC",
                   "PTGES","LRP2","HNF1B",
                   "ARGFX","TDGF1","SLC34A2",
                   "PDGFA","P4HA1")
for(i in selectedGenes){
  pdf(paste('./ESpseudotime',i,'.pdf', sep=''),width = 3.5, height = 2)
  print(ggplot(ESCscaledData,
               aes_string(x="pseudotime",y=i)) + geom_point(shape=19, size=1.5, color='#feb24c') +
          geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=i), color='#54278f')+
          ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
          geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme)
  dev.off()
}

mytheme <-  theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 1,
                                           colour = "black"),
                  axis.title.x =element_blank(),
                  axis.text.x = element_text(size = 15,
                                             family = "sans",
                                             color = "black",
                                          
                                             vjust = 0,
                                             hjust = 0),
                  axis.text.y = element_text(size = 15,
                                             family = "sans",
                                             color = "black",
                                            
                                             vjust = 0,
                                             hjust = 1),
                  axis.title.y=element_blank(),
                  legend.text = element_text(size=15,
                                             family = "sans",
                                             color = "black"),
                  legend.title = element_text(size=15,
                                              family = "sans",
                                              color = "black"),
                  legend.background = element_blank(),
                  legend.key=element_blank(),
                  plot.title=element_text(family="sans",size=15,color="black",hjust=0.5,lineheight=0.5,vjust=0.5))
selectedGenes <- c("GATA2","PDGFA","PTGES",
                   "MYC","SLC34A2","ARGFX")
x <- c(1,2,3,1,2,3)
y <- c(1,1,1,2,2,2)
pdf(file = './combine_allscatters.pdf', width = 9, height = 6)
grid.newpage()  ###新建图表版面
pushViewport(viewport(layout = grid.layout(3,2))) ####将版面分面
vplayout <- function(x,y){
  viewport(layout.pos.row = x, layout.pos.col = y)
}
for(i in 1:length(selectedGenes)){
print(ggplot(ESCscaledData,
                   aes_string(x="pseudotime",y=selectedGenes[i])) + geom_point(shape=19, size=1.5, color='#feb24c') +
              geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[i]), color='#54278f')+
              ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[i])+
              geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme , vp = vplayout(x[i],y[i]))
}
dev.off()

p1 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[1])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[1]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[1])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p2 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[2])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[2]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[2])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p3 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[3])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[3]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[3])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p4 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[4])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[4]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[4])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p5 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[5])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[5]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[5])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
p6 <- ggplot(ESCscaledData,
       aes_string(x="pseudotime",y=selectedGenes[6])) + geom_point(shape=19, size=1.5, color='#feb24c') +
  geom_point(data=transitionPoint_plot,shape=19, aes_string(x="pseudotime",y=selectedGenes[6]), color='#54278f')+
  ylab("scaled data") + xlab("pseudo time") + ggtitle(selectedGenes[6])+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+mytheme 
pdf(file = './combine_allscatters.pdf', width = 6, height = 9)
ggpubr::ggarrange(
  p1,p5,p3,p4,p2,p6,
  ncol = 2, nrow = 3,align = 'v')
dev.off()  
  
#plot the path
basin[which(basin==FALSE)] <- 0
basin[which(basin=='transition path\nfrom 3 to 8')] <- 1
basin[which(basin=='transitionPoint')] <- 2
basin <- as.numeric(basin)
pdf(file = './ESPCpathplot.pdf')
scatter3D(x=phate.ESCobject$embedding[,1],
          y=phate.ESCobject$embedding[,2], 
          z=phate.ESCobject$embedding[,3],
          colvar = basin, 
          col = c(alpha("#737373",0.1),
            alpha("#feb24c",1),
            alpha('#54278f',1)),
          pch = 19, 
          cex = 1,
          theta = 0,
          phi = 125)
dev.off()

save.image('./ESCDevelopmentPath.RData')

```
