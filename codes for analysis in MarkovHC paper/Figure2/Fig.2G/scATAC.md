```R
setwd('/data02/zywang/MarkovHC/supplementaryFigures/scRNAandscATAC/')
set.seed(1234)
options(stringsAsFactors = F)
library(flowCore)
library(Rtsne)
library(ggplot2)
library(MarkovHC)
library(RColorBrewer)
library(pheatmap)
library(stringr)
library(phateR)
library(plot3D)
library(hdf5r)
library(stringi)
library(mclust)
library(aricode)
library(Seurat)
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


```R
lsi_data <- readRDS('./lsi_data.Rds')
peak_counts <- readRDS('./peak_counts.Rds')
```


```R
#realLabels are the real labels of each sample.
#comparedMethods is a character vector method names.
realLabels <- readRDS('./realLabels.Rds')
evaluation_dataFrame <- as.data.frame(matrix(0,9,2))
colnames(evaluation_dataFrame) <- c('ARI', 'NMI')
rownames(evaluation_dataFrame) <- c('MarkovHC','Seurat','SIMLR','SC3','kmeans','HC','hdbscan','specc', 'mclust')
```

# SIMLR


```R
dim(lsi_data)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>8000</li><li>39</li></ol>



# SIMLR always collapses when analyses large datasets, so we sampled 8k cells from this dataset.


```R
SIMLRObject = SIMLR(X =  lsi_data%>%t(), c = 29)
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
    Performing t-SNE.
    Epoch: Iteration # 100  error is:  0.8216889 
    Epoch: Iteration # 200  error is:  0.6695331 
    Epoch: Iteration # 300  error is:  0.6125788 
    Epoch: Iteration # 400  error is:  0.5827744 
    Epoch: Iteration # 500  error is:  0.5643749 
    Epoch: Iteration # 600  error is:  0.55184 
    Epoch: Iteration # 700  error is:  0.5426433 
    Epoch: Iteration # 800  error is:  0.5356395 
    Epoch: Iteration # 900  error is:  0.5304242 
    Epoch: Iteration # 1000  error is:  0.5266997 
    Performing Kmeans.
    Performing t-SNE.
    Epoch: Iteration # 100  error is:  14.09709 
    Epoch: Iteration # 200  error is:  0.951253 
    Epoch: Iteration # 300  error is:  0.8167509 
    Epoch: Iteration # 400  error is:  0.7672013 
    Epoch: Iteration # 500  error is:  0.7419637 
    Epoch: Iteration # 600  error is:  0.7271148 
    Epoch: Iteration # 700  error is:  0.7176807 
    Epoch: Iteration # 800  error is:  0.7113918 
    Epoch: Iteration # 900  error is:  0.7074404 
    Epoch: Iteration # 1000  error is:  0.7047975 


# sc3


```R
sce <- SingleCellExperiment(
assays = list(
    counts = as.matrix(peak_counts),
    logcounts = as.matrix(peak_counts)
    )
)
rowData(sce)$feature_symbol <- rownames(peak_counts)
sce <- sc3(sce, ks = 29, biology = FALSE, gene_filter = FALSE)
```

    Setting SC3 parameters...
    
    Your dataset contains more than 2000 cells. Adjusting the nstart parameter of kmeans to 50 for faster performance...
    
    Defining training cells for SVM using 5000 random cells...
    
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
kmeans_results <- kmeans(lsi_data, centers=29)
```


```R
evaluation_dataFrame$ARI[5] <- adjustedRandIndex(realLabels, as.character(kmeans_results$cluster))
evaluation_dataFrame$NMI[5] <- NMI(realLabels, as.character(kmeans_results$cluster))
```

# hierarchical average


```R
hresult_average <- hclust(dist(lsi_data),method = 'average')
hresult_average <- cutree(hresult_average, k=29)
```


```R
evaluation_dataFrame$ARI[6] <- adjustedRandIndex(realLabels, as.character(hresult_average))
evaluation_dataFrame$NMI[6] <- NMI(realLabels, as.character(hresult_average))
```

# hdbscan


```R
hdbscan_res <- hdbscan(lsi_data, minPts=10)
hdbscan_res <- hdbscan_res$cluster
```


```R
evaluation_dataFrame$ARI[7] <- adjustedRandIndex(realLabels, as.character(hdbscan_res))
evaluation_dataFrame$NMI[7] <- NMI(realLabels, as.character(hdbscan_res))
```

# specc


```R
sp_result <- kernlab::specc(lsi_data, centers=29)
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
EM_res <- mclust::Mclust( lsi_data )
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
	<tr><th scope=row>MarkovHC</th><td>0.000000000</td><td>0.00000000</td></tr>
	<tr><th scope=row>Seurat</th><td>0.000000000</td><td>0.00000000</td></tr>
	<tr><th scope=row>SIMLR</th><td>0.347310217</td><td>0.53683187</td></tr>
	<tr><th scope=row>SC3</th><td>0.187245808</td><td>0.41055077</td></tr>
	<tr><th scope=row>kmeans</th><td>0.269591655</td><td>0.49775413</td></tr>
	<tr><th scope=row>HC</th><td>0.003748962</td><td>0.01839170</td></tr>
	<tr><th scope=row>hdbscan</th><td>0.114878622</td><td>0.09330571</td></tr>
	<tr><th scope=row>specc</th><td>0.245372547</td><td>0.49811439</td></tr>
	<tr><th scope=row>mclust</th><td>0.622057981</td><td>0.58164994</td></tr>
</tbody>
</table>




```R
saveRDS(evaluation_dataFrame, './evaluation_dataFrame_scATAC.RDs')
```


```R
save(
    SIMLRObject,
    sc_labels,
    kmeans_results,
    hresult_average,
    hdbscan_res,
    sp_result,
    EM_res,
    file = './scATAC_otherMethods.RData')
```


```R
save.image('./scATAC.RData')
```


```R

```
