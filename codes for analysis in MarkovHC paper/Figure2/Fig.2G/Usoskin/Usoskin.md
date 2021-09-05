```R
library(SingleCellExperiment)
library(SC3)
library(scater)
library(Seurat)
library(MarkovHC)
library(ggplot2)
library(EMCluster)
library(cluster)
library(dplyr)
library(mclust)
library(reshape2)
library(dbscan)
library(SIMLR)
library(aricode)
library(Hmisc)
library(clusterProfiler)
library(stringr)
options(repr.plot.width=5, repr.plot.height=5)
setwd('/data02/zywang/MarkovHC/Figure3/')
```


```R
#Figures
mytheme <-  theme(panel.grid.major =element_blank(),
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
                                             vjust = 1,
                                             hjust = 1,
                                            angle=45),
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
                  legend.key=element_blank(),
                  plot.title=element_text(family="sans",size=15,color="black",
                                          face="bold",hjust=0.5,lineheight=0.5,vjust=0.5))

notheme <- mytheme+NoLegend()
```

# load data


```R
load('./Test_Usoskin.RData')
```


```R
rownames(Test_4_Usoskin[[1]]) <- paste('gene',1:nrow(Test_4_Usoskin[[1]]),sep='')
colnames(Test_4_Usoskin[[1]]) <- paste('cell',1:ncol(Test_4_Usoskin[[1]]),sep='')
```


```R
input_matrix <- Test_4_Usoskin[[1]]
```

# preprocessing


```R
SeuratObject <- CreateSeuratObject(counts=input_matrix, 
                                   project = "Usoskin",
                                   min.cells = 0, 
                                   min.features = 0)
SeuratObject <- SetAssayData(object = SeuratObject, 
                             slot = "scale.data", 
                             new.data = input_matrix)
```


```R
SeuratObject@meta.data$label <- Test_4_Usoskin[[3]][,1]
```


```R
SeuratObject <- RunPCA(SeuratObject, 
                       npcs = 100,
                       features = rownames(SeuratObject), 
                       verbose=FALSE)
```


```R
ElbowPlot(SeuratObject, ndims = 100)
```


![png](output_10_0.png)


# PC selection


```R
PC_selection(SeuratObject)
```

    [1] 4



![png](output_12_1.png)



```R
SeuratObject <- RunUMAP(object = SeuratObject, dims=1:4)
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    19:29:11 UMAP embedding parameters a = 0.9922 b = 1.112
    
    19:29:11 Read 622 rows and found 4 numeric columns
    
    19:29:11 Using Annoy for neighbor search, n_neighbors = 30
    
    19:29:11 Building Annoy index with metric = cosine, n_trees = 50
    
    0%   10   20   30   40   50   60   70   80   90   100%
    
    [----|----|----|----|----|----|----|----|----|----|
    
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    *
    
    |
    
    19:29:11 Writing NN index file to temp file /tmp/Rtmp4Iqn4G/file46d5b365467a
    
    19:29:11 Searching Annoy index using 1 thread, search_k = 3000
    
    19:29:11 Annoy recall = 100%
    
    19:29:12 Commencing smooth kNN distance calibration using 1 thread
    
    19:29:13 Initializing from normalized Laplacian + noise
    
    19:29:13 Commencing optimization for 500 epochs, with 20762 positive edges
    
    19:29:15 Optimization finished
    



```R
SeuratObject <- FindNeighbors(object = SeuratObject,
                              k.param = 20,
                              compute.SNN = TRUE,
                              prune.SNN = 0,
                              reduction = "pca", 
                              dims = 1:4,
                              force.recalc = TRUE)
```

    Computing nearest neighbor graph
    
    Computing SNN
    



```R
DimPlot(SeuratObject, reduction = "umap", group.by = 'label')
```


![png](output_15_0.png)



```R
#realLabels are the real labels of each sample.
#comparedMethods is a character vector method names.
realLabels=SeuratObject$label
comparedMethods=c('MarkovHC','Seurat','SIMLR','SC3','kmeans','HC','hdbscan','specc', 'mclust')
```


```R
evaluation_dataFrame <- as.data.frame(matrix(0, nrow = length(comparedMethods), ncol = 2))
rownames(evaluation_dataFrame) <- comparedMethods
colnames(evaluation_dataFrame) <- c('ARI', 'NMI')
```

# run MarkovHC


```R
MarkovHC_object <- MarkovHC(MarkovHC_input = SeuratObject,
                            dobasecluster = FALSE,
                            SNNslot = 'RNA_snn', 
                            KNNslot = 'RNA_nn',
                            cutpoint = 0.001,
                            verbose = FALSE)
```

    [1] "The input is a Seurat object."


# level selection


```R
energyGap_selection(MarkovObject=MarkovHC_object, m=3)
```

    [1] "levels with possible biological meaning:"
    0.1% 0.2% 0.3% 0.6%  50% 
      12   21   29   37   50 
    [1] "the level may with an optimal cluster number is among:"
    [1] "levels:from 46 to 50"



![png](output_21_1.png)



```R
internal_measures <- IMI_selection(MarkovObject=MarkovHC_object,
                                   prune=TRUE,
                                   weed=10)
```


```R
head(internal_measures, n=10)
```


<table>
<caption>A data.frame: 10 × 6</caption>
<thead>
	<tr><th></th><th scope=col>Name</th><th scope=col>Score</th><th scope=col>connectivity</th><th scope=col>silhouette</th><th scope=col>dunn</th><th scope=col>C_cut_gap</th></tr>
	<tr><th></th><th scope=col>&lt;int&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>50</th><td>50</td><td>0.001419618</td><td>  73.68416</td><td>0.43943475</td><td>0.004853957</td><td>13.78386652</td></tr>
	<tr><th scope=row>48</th><td>48</td><td>0.005912611</td><td>  58.75982</td><td>0.47433161</td><td>0.004853957</td><td> 1.80685575</td></tr>
	<tr><th scope=row>2</th><td> 2</td><td>0.007265416</td><td> 203.62712</td><td>0.02955078</td><td>0.007489736</td><td> 0.86584422</td></tr>
	<tr><th scope=row>3</th><td> 3</td><td>0.023754505</td><td>2638.34159</td><td>0.32126354</td><td>0.005046630</td><td> 0.36131894</td></tr>
	<tr><th scope=row>49</th><td>49</td><td>0.038748886</td><td>  45.56157</td><td>0.36147756</td><td>0.004853957</td><td> 2.18612664</td></tr>
	<tr><th scope=row>42</th><td>42</td><td>0.077053734</td><td>  71.06495</td><td>0.37228058</td><td>0.004759681</td><td> 0.14582985</td></tr>
	<tr><th scope=row>46</th><td>46</td><td>0.102879427</td><td>  93.85368</td><td>0.36626549</td><td>0.004759681</td><td> 0.43445755</td></tr>
	<tr><th scope=row>35</th><td>35</td><td>0.138506454</td><td>  15.18164</td><td>0.35366662</td><td>0.004759681</td><td> 0.09307756</td></tr>
	<tr><th scope=row>7</th><td> 7</td><td>0.165459090</td><td>  16.33283</td><td>0.36537797</td><td>0.003862940</td><td> 0.16467321</td></tr>
	<tr><th scope=row>31</th><td>31</td><td>0.201634807</td><td>  81.71917</td><td>0.35480677</td><td>0.004947639</td><td> 0.04453450</td></tr>
</tbody>
</table>




```R
MarkovHCLabels <-  fetchLabels(MarkovObject=MarkovHC_object,
                               MarkovLevels=3:length(MarkovHC_object$hierarchicalStructure),
                               prune = TRUE, weed = 10)
```


```R
length(unique( MarkovHCLabels$lv48))
```


4



```R
MarkovHCLabels <-  MarkovHCLabels$lv48
```


```R
evaluation_dataFrame$ARI[1] <- adjustedRandIndex(realLabels, MarkovHCLabels)
evaluation_dataFrame$NMI[1] <- NMI(realLabels, MarkovHCLabels)
```

# Seurat


```R
SeuratObject <- FindClusters(SeuratObject)
```

    Modularity Optimizer version 1.3.0 by Ludo Waltman and Nees Jan van Eck
    
    Number of nodes: 622
    Number of edges: 20869
    
    Running Louvain algorithm...
    Maximum modularity in 10 random starts: 0.7944
    Number of communities: 8
    Elapsed time: 0 seconds



```R
evaluation_dataFrame$ARI[2] <- adjustedRandIndex(realLabels, as.character(SeuratObject@meta.data$seurat_clusters))
evaluation_dataFrame$NMI[2] <- NMI(realLabels, as.character(SeuratObject@meta.data$seurat_clusters))   
```

# SIMLR


```R
SIMLRObject = SIMLR(X =  Embeddings(object = SeuratObject, reduction = "pca")[,1:4]%>%t(), 
                    c = 4)
evaluation_dataFrame$ARI[3] <- adjustedRandIndex(realLabels, as.character(SIMLRObject$y$cluster))
evaluation_dataFrame$NMI[3] <- NMI(realLabels, as.character(SIMLRObject$y$cluster))    
```

    Computing the multiple Kernels.
    Performing network diffiusion.
    Iteration:  1 
    Iteration:  2 
    Iteration:  3 
    Iteration:  4 
    Iteration:  5 
    Iteration:  6 
    Iteration:  7 
    Iteration:  8 
    Iteration:  9 
    Iteration:  10 
    Iteration:  11 
    Iteration:  12 
    Iteration:  13 
    Iteration:  14 
    Iteration:  15 
    Iteration:  16 
    Iteration:  17 
    Iteration:  18 
    Performing t-SNE.
    Epoch: Iteration # 100  error is:  0.1897101 
    Epoch: Iteration # 200  error is:  0.1266257 
    Epoch: Iteration # 300  error is:  0.1207027 
    Epoch: Iteration # 400  error is:  0.1176287 
    Epoch: Iteration # 500  error is:  0.1156548 
    Epoch: Iteration # 600  error is:  0.114244 
    Epoch: Iteration # 700  error is:  0.1131753 
    Epoch: Iteration # 800  error is:  0.112319 
    Epoch: Iteration # 900  error is:  0.1116131 
    Epoch: Iteration # 1000  error is:  0.1110174 
    Performing Kmeans.
    Performing t-SNE.
    Epoch: Iteration # 100  error is:  10.74308 
    Epoch: Iteration # 200  error is:  0.1879954 
    Epoch: Iteration # 300  error is:  0.1538846 
    Epoch: Iteration # 400  error is:  0.1473543 
    Epoch: Iteration # 500  error is:  0.1441514 
    Epoch: Iteration # 600  error is:  0.1421812 
    Epoch: Iteration # 700  error is:  0.1408088 
    Epoch: Iteration # 800  error is:  0.1397686 
    Epoch: Iteration # 900  error is:  0.1389501 
    Epoch: Iteration # 1000  error is:  0.1382885 


# sc3


```R
sce <- SingleCellExperiment(
assays = list(
    counts = as.matrix(GetAssayData(object = SeuratObject, slot = "counts")),
    logcounts = as.matrix(GetAssayData(object = SeuratObject, slot = "counts"))
    )
)
rowData(sce)$feature_symbol <- rownames(GetAssayData(object = SeuratObject, slot = "counts"))
sce <- sc3(sce, ks = 4, biology = FALSE)
```

    Setting SC3 parameters...
    
    Warning message:
    “'isSpike' is deprecated.
    See help("Deprecated")”
    Calculating distances between the cells...
    
    Performing transformations and calculating eigenvectors...
    
    Performing k-means clustering...
    


    


    Calculating consensus matrix...
    



```R
sc_labels <- as.character(sce@colData[,1])
sc_labels[which(is.na(sc_labels))] <- "0"
evaluation_dataFrame$ARI[4] <- adjustedRandIndex(realLabels, sc_labels)
evaluation_dataFrame$NMI[4] <- NMI(realLabels, sc_labels)  
```

# kmeans


```R
kmeans_results <- kmeans(Embeddings(object = SeuratObject, reduction = "pca")[,1:4], centers=4)
```


```R
evaluation_dataFrame$ARI[5] <- adjustedRandIndex(realLabels, as.character(kmeans_results$cluster))
evaluation_dataFrame$NMI[5] <- NMI(realLabels, as.character(kmeans_results$cluster))
```

# hierarchical average


```R
hresult_average <- hclust(dist(Embeddings(object = SeuratObject, reduction = "pca")[,1:4]),method = 'average')
hresult_average <- cutree(hresult_average, k=4)
```


```R
evaluation_dataFrame$ARI[6] <- adjustedRandIndex(realLabels, as.character(hresult_average))
evaluation_dataFrame$NMI[6] <- NMI(realLabels, as.character(hresult_average))
```

# hdbscan


```R
hdbscan_res <- hdbscan(Embeddings(object = SeuratObject, reduction = "pca")[,1:4], minPts=10)
hdbscan_res <- hdbscan_res$cluster
```


```R
evaluation_dataFrame$ARI[7] <- adjustedRandIndex(realLabels, as.character(hdbscan_res))
evaluation_dataFrame$NMI[7] <- NMI(realLabels, as.character(hdbscan_res))
```

# specc


```R
sp_result <- kernlab::specc(Embeddings(object = SeuratObject, reduction = "pca")[,1:4], centers=4)
```


```R
sp_result <- sp_result@.Data
```


```R
evaluation_dataFrame$ARI[8] <- adjustedRandIndex(realLabels, as.character(sp_result))
evaluation_dataFrame$NMI[8] <- NMI(realLabels, as.character(sp_result))
```

# mclust


```R
EM_res <- mclust::Mclust( Embeddings(object = SeuratObject, reduction = "pca")[,1:4] )
```


```R
evaluation_dataFrame$ARI[9] <- adjustedRandIndex(realLabels, as.character(EM_res$classification))
evaluation_dataFrame$NMI[9] <- NMI(realLabels, as.character(EM_res$classification))
```


```R
evaluation_dataFrame
```


<table>
<caption>A data.frame: 9 × 2</caption>
<thead>
	<tr><th></th><th scope=col>ARI</th><th scope=col>NMI</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>MarkovHC</th><td>0.6756813</td><td>0.6524085</td></tr>
	<tr><th scope=row>Seurat</th><td>0.3669714</td><td>0.4445319</td></tr>
	<tr><th scope=row>SIMLR</th><td>0.7535641</td><td>0.6521422</td></tr>
	<tr><th scope=row>SC3</th><td>0.8848992</td><td>0.8667264</td></tr>
	<tr><th scope=row>kmeans</th><td>0.4032043</td><td>0.4262882</td></tr>
	<tr><th scope=row>HC</th><td>0.2682046</td><td>0.3176759</td></tr>
	<tr><th scope=row>hdbscan</th><td>0.5344517</td><td>0.5551469</td></tr>
	<tr><th scope=row>specc</th><td>0.7805418</td><td>0.7031934</td></tr>
	<tr><th scope=row>mclust</th><td>0.4965239</td><td>0.5536991</td></tr>
</tbody>
</table>



# Figures


```R
allColors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3","#ff7f00","#ffff33","#a65628","#f781bf","#999999","#8dd3c7","#ffffb3","#bebada","#fb8072","#80b1d3","#fdb462","#b3de69","#fccde5","#a6cee3","#1f78b4","#b2df8a",
"#33a02c","#fb9a99","#e31a1c","#fdbf6f","#cab2d6","#fbb4ae","#b3cde3","#ccebc5","#decbe4","#fed9a6","#ffffcc","#e5d8bd","#fddaec","#8c96c6","#8c6bb1","#88419d","#810f7c","#4d004b","#74c476","#41ab5d",
"#238b45","#006d2c","#00441b","#fe9929","#ec7014","#cc4c02","#993404","#662506","#df65b0","#e7298a","#ce1256","#980043","#67001f")
```


```R
SeuratObject@meta.data$MarkovHC <- MarkovHCLabels
SeuratObject@meta.data$SIMLR <- as.character(SIMLRObject$y$cluster)
SeuratObject@meta.data$SC3 <- sc_labels
SeuratObject@meta.data$kmeans <- as.character(kmeans_results$cluster)
SeuratObject@meta.data$HC <- as.character(hresult_average)
SeuratObject@meta.data$hdbscan <- as.character(hdbscan_res)
SeuratObject@meta.data$specc <- as.character(sp_result)
SeuratObject@meta.data$mclust <- as.character(as.character(EM_res$classification))
```


```R
colorSet = function(seuratObject=NULL,
                    colorVector=NULL,
                    method=NULL){
    seuratObject@meta.data[,method] <- as.character(seuratObject@meta.data[,method])
    label2label <- as.data.frame(unique(seuratObject@meta.data[,method]),
                                 stringsAsFactors = FALSE)
    label2label$V2 <- label2label[,1]
    for(i in label2label[,1]){
        temp <- subset(seuratObject@meta.data, seuratObject@meta.data[,method]==i)
        tempLabel <- temp$label
        tempLabel_feq <- table(tempLabel)
        label2label[which(label2label[,1]==i),2] <- as.numeric(names(tempLabel_feq)[tempLabel_feq == max(tempLabel_feq)])[1]
    }
    colors <- colorVector[as.numeric(label2label[,2])]
    colors_fre <- table(colors)
    repeatcolors <- names(colors_fre)[colors_fre >1] 
    colors[which(colors%in%repeatcolors)] <- sample(allColors,length(which(colors%in%repeatcolors)))
    names(colors) <- label2label[,1]
    return(colors)
}
```


```R
colorVector <-  c('#e41a1c','#377eb8','#4daf4a','#984ea3')
```


```R
for(i in c('MarkovHC','seurat_clusters','SIMLR','SC3','kmeans','HC','hdbscan','specc','mclust')){
    colorVector.temp <- colorSet(seuratObject=SeuratObject,
                                 colorVector=colorVector,
                                 method=i)
    assign(paste(i,'_plot_Usoskin',sep=''), value = DimPlot(SeuratObject, group.by=i, cols=colorVector.temp, pt.size=2)+notheme)
}
```


```R
names(colorVector) <- 1:length(unique(SeuratObject@meta.data$label))
```


```R
groundTruth_plot_Usoskin <- DimPlot(SeuratObject, group.by="label", cols=colorVector, pt.size=2)+notheme
```


```R
save(
    groundTruth_plot_Usoskin,
    MarkovHC_plot_Usoskin,
    seurat_clusters_plot_Usoskin,
    SIMLR_plot_Usoskin,
    SC3_plot_Usoskin,
    kmeans_plot_Usoskin,
    HC_plot_Usoskin,
    hdbscan_plot_Usoskin,
    specc_plot_Usoskin,
    mclust_plot_Usoskin,
    file = './Usoskin_plot.RData')
```


```R
saveRDS(evaluation_dataFrame, './evaluation_dataFrame_Usoskin.RDs')
```


```R
save.image('./Usoskin.RData')
```


```R

```
