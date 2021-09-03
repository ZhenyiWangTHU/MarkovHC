```R
setwd('/data02/zywang/MarkovHC/Figure6/')
library(MarkovHC)
library(ggplot2)
library(Seurat)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MASS)
library(stringr)
library(pheatmap)
library(phateR)
library(plot3D)
library(monocle)
```


```R
load('./dt.annotation.v20190230.Rdata')
```


```R
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
                  legend.key=element_blank(),
                  plot.title=element_text(family="sans",size=15,color="black",
                                          face="bold",hjust=0.5,lineheight=0.5,vjust=0.5))
```


```R
dt
```


    An old seurat object
     22910 genes across 32332 samples



```R
dt_update <- UpdateSeuratObject(dt)
```

    Updating from v2.X to v3.X
    
    Validating object structure
    
    Updating object slots
    
    Ensuring keys are in the proper strucutre
    
    Ensuring feature names don't have underscores or pipes
    
    Object representation is consistent with the most current Seurat version
    



```R
ElbowPlot(dt_update, ndims = 50)
```


![png](output_5_0.png)



```R
DimPlot(dt_update, reduction = "tsne", label = TRUE)
```

    Warning message:
    “Using `as.character()` on a quosure is deprecated as of rlang 0.3.0.
    Please use `as_label()` or `as_name()` instead.
    [90mThis warning is displayed once per session.[39m”



![png](output_6_1.png)



```R
dt_update <- RunUMAP(object = dt_update, dims=1:20, n.neighbors = 50, min.dist = 1, seed.use =1)
```

    Warning message:
    “The default method for RunUMAP has changed from calling Python UMAP via reticulate to the R-native UWOT using the cosine metric
    To use Python UMAP via reticulate, set umap.method to 'umap-learn' and metric to 'correlation'
    This message will be shown once per session”
    15:54:25 UMAP embedding parameters a = 0.115 b = 1.929
    
    15:54:25 Read 32332 rows and found 20 numeric columns
    
    15:54:25 Using Annoy for neighbor search, n_neighbors = 50
    
    15:54:25 Building Annoy index with metric = cosine, n_trees = 50
    
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
    
    15:54:30 Writing NN index file to temp file /tmp/RtmpVLeiSf/file252b521e46225
    
    15:54:30 Searching Annoy index using 1 thread, search_k = 5000
    
    15:54:47 Annoy recall = 100%
    
    15:54:47 Commencing smooth kNN distance calibration using 1 thread
    
    15:54:50 Initializing from normalized Laplacian + noise
    
    15:54:53 Commencing optimization for 200 epochs, with 2351666 positive edges
    
    15:55:38 Optimization finished
    



```R
DimPlot(dt_update, reduction = "umap", label = TRUE)
```


![png](output_8_0.png)



```R
unique(dt_update@meta.data$batch)
#MSC and Cancer cell of EGC
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>'NAG1'</li><li>'NAG2'</li><li>'NAG3'</li><li>'CAG1'</li><li>'CAG2'</li><li>'CAG3'</li><li>'IMW1'</li><li>'IMW2'</li><li>'IMS1'</li><li>'IMS2'</li><li>'IMS3'</li><li>'IMS4'</li><li>'EGC'</li></ol>




```R
dt_sub <- subset(dt_update, subset = batch == 'EGC')
```


```R
unique(Idents(dt_sub))
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>MSC</li><li>Enteroendocrine</li><li>Cancer cell</li><li>PMC</li><li>PC</li><li>Enterocyte</li><li>T cell</li><li>Goblet cell</li><li>B cell</li><li>Macrophage</li><li>EC</li><li>GMC</li><li>Chief cell</li><li>Mast cell</li><li>Fibroblast</li></ol>

<details>
	<summary style=display:list-item;cursor:pointer>
		<strong>Levels</strong>:
	</summary>
	<style>
	.list-inline {list-style: none; margin:0; padding: 0}
	.list-inline>li {display: inline-block}
	.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
	</style>
	<ol class=list-inline><li>'PMC'</li><li>'GMC'</li><li>'Chief cell'</li><li>'B cell'</li><li>'Cancer cell'</li><li>'EC'</li><li>'Enterocyte'</li><li>'Enteroendocrine'</li><li>'Fibroblast'</li><li>'Goblet cell'</li><li>'Macrophage'</li><li>'Mast cell'</li><li>'MSC'</li><li>'PC'</li><li>'T cell'</li></ol>
</details>



```R
dt_sub <- subset(dt_sub, idents=c('PMC','Cancer cell', 'Enterocyte',
                                  'MSC', 'Chief cell',
                                  'PC', 'GMC'))
```


```R
dt_sub <- subset(dt_sub, idents=c('Cancer cell','MSC'))
```


```R
dt_sub
```


    An object of class Seurat 
    22910 features across 1526 samples within 1 assay 
    Active assay: RNA (22910 features, 1258 variable features)
     3 dimensional reductions calculated: pca, tsne, umap



```R
dt_sub <- NormalizeData(dt_sub, normalization.method = "LogNormalize", scale.factor = 10000)
```


```R
dt_sub <- FindVariableFeatures(dt_sub, selection.method = "vst", nfeatures = 3000)
# Identify the 10 most highly variable genes
dt_sub10 <- head(VariableFeatures(dt_sub), 10)
# plot variable features with and without labels
dt_sub1 <- VariableFeaturePlot(dt_sub)
dt_sub2 <- LabelPoints(plot = dt_sub1, points = dt_sub10, repel = TRUE)
dt_sub2
```

    When using repel, set xnudge and ynudge to 0 for optimal results
    
    Warning message:
    “Transformation introduced infinite values in continuous x-axis”



![png](output_16_1.png)



```R
dt_sub <- ScaleData(dt_sub, features = rownames(dt_sub) ,vars.to.regress =  NULL)
```

    Centering and scaling data matrix
    



```R
dt_sub <- RunPCA(dt_sub,  
                 npcs = 1000, 
                 features = VariableFeatures(object = dt_sub), 
                 verbose=FALSE)
```

    Warning message in irlba(A = t(x = object), nv = npcs, ...):
    “You're computing too large a percentage of total singular values, use a standard svd instead.”



```R
ElbowPlot(dt_sub,ndims = 1000)
```


![png](output_19_0.png)


# PC selection


```R
PC_selection(dt_sub)
```

    [1] 8



![png](output_21_1.png)



```R
dt_sub <- RunUMAP(object = dt_sub, 
                  reduction = "pca", 
                  dims=1:8,
                  n.neighbors=60L, 
                  min.dist=0.6, 
                  seed.use=100L, 
                  n.components=2,
                  umap.method = 'umap-learn', metric = 'correlation')
```


```R
DimPlot(dt_sub, reduction = "umap", label = TRUE,label.size=10,dims = c(1, 2))+NoLegend()
```


![png](output_23_0.png)



```R
layout <- as.data.frame(Embeddings(object = dt_sub, reduction = "umap"))
layout$idents <- Idents(dt_sub)
```


```R
head(layout)
```


<table>
<caption>A data.frame: 6 × 3</caption>
<thead>
	<tr><th></th><th scope=col>UMAP_1</th><th scope=col>UMAP_2</th><th scope=col>idents</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>AAACCTGAGAAACGAG-14</th><td> 3.6527257</td><td>-0.5297089</td><td>MSC        </td></tr>
	<tr><th scope=row>AAACCTGAGGGCACTA-14</th><td> 1.7683476</td><td> 2.5783837</td><td>MSC        </td></tr>
	<tr><th scope=row>AAACCTGCAATGGTCT-14</th><td>-7.3420510</td><td> 4.9964356</td><td>MSC        </td></tr>
	<tr><th scope=row>AAACCTGCATGGGAAC-14</th><td> 0.9668099</td><td>-2.2345328</td><td>MSC        </td></tr>
	<tr><th scope=row>AAACCTGCATGGTTGT-14</th><td> 0.6070247</td><td>-1.2994579</td><td>MSC        </td></tr>
	<tr><th scope=row>AAACCTGGTTAGTGGG-14</th><td>-7.3071737</td><td>-1.6284339</td><td>Cancer cell</td></tr>
</tbody>
</table>




```R
pdf(file = './groundtruth_nolegend.pdf', width = 3.5, height = 3.5)
ggplot(data=layout, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=idents), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
  scale_fill_manual(
    values = c( "Cancer cell"=alpha("#cb181d",0.7),
                "MSC"=alpha("#08519c",0.7)),
    breaks = c( "Cancer cell",
                 "MSC"))
dev.off()
```


<strong>png:</strong> 2



```R
table(Idents(dt_sub))
```


    
    Cancer cell         MSC 
            695         831 



```R
dt_sub <- FindNeighbors(object = dt_sub,
                        k.param = 60,
                        compute.SNN = TRUE,
                        prune.SNN = 0,
                        reduction = "pca", 
                        dims = 1:8,
                        force.recalc = TRUE)
```

    Computing nearest neighbor graph
    
    Computing SNN
    


# run MarkovHC


```R
MarkovHC_dt <- MarkovHC(dt_sub,
                        dobasecluster = TRUE,
                        cutpoint = 0.001,
                        verbose = FALSE)
```

    [1] "The input is a Seurat object."


# level selection


```R
energyGap_selection(MarkovObject=MarkovHC_dt, m=3)
```

    [1] "levels with possible biological meaning:"
    0.5% 0.8% 1.8%  50% 
       6   10   14   20 
    [1] "the level may with an optimal cluster number is among:"
    [1] "levels:from 9 to 20"



![png](output_32_1.png)



```R
internal_measures <- IMI_selection(MarkovObject=MarkovHC_dt,
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
	<tr><th scope=row>18</th><td>18</td><td>0.01120142</td><td> 8.350000</td><td>-0.186499137</td><td>0.007889161</td><td>0.30284685</td></tr>
	<tr><th scope=row>20</th><td>20</td><td>0.02212964</td><td> 3.783333</td><td> 0.004008972</td><td>0.007889161</td><td>0.89155823</td></tr>
	<tr><th scope=row>21</th><td>21</td><td>0.03642169</td><td>10.750000</td><td>-0.037341993</td><td>0.010612066</td><td>0.01452397</td></tr>
	<tr><th scope=row>17</th><td>17</td><td>0.06994058</td><td> 2.866667</td><td>-0.148025095</td><td>0.007656491</td><td>0.21517820</td></tr>
	<tr><th scope=row>14</th><td>14</td><td>0.25000000</td><td> 2.500000</td><td>-0.199536861</td><td>0.007656491</td><td>0.21729775</td></tr>
	<tr><th scope=row>2</th><td> 2</td><td>0.48768868</td><td> 7.083333</td><td>-0.195952812</td><td>0.004514750</td><td>0.09869341</td></tr>
	<tr><th scope=row>16</th><td>16</td><td>0.48768868</td><td> 2.033333</td><td>-0.157311424</td><td>0.007656491</td><td>0.06371887</td></tr>
	<tr><th scope=row>19</th><td>19</td><td>0.48768868</td><td> 1.916667</td><td>-0.174001146</td><td>0.007889161</td><td>0.06504726</td></tr>
	<tr><th scope=row>3</th><td> 3</td><td>0.86443890</td><td> 2.000000</td><td>-0.151512100</td><td>0.004514750</td><td>0.03712145</td></tr>
	<tr><th scope=row>13</th><td>13</td><td>1.00000000</td><td> 5.566667</td><td>-0.242616175</td><td>0.005256388</td><td>0.02833771</td></tr>
</tbody>
</table>




```R
MarkovHCLabels <-  fetchLabels(MarkovObject=MarkovHC_dt,
                               MarkovLevels=1:length(MarkovHC_dt$hierarchicalStructure),
                               prune = TRUE, weed = 10)
```


```R
length(unique(MarkovHCLabels$lv18))
```


5


# Figure


```R
DimPlot(dt_sub, reduction = "umap", label = TRUE,label.size=10,dims = c(1, 2))+NoLegend()
```


![png](output_38_0.png)


# lv18

# transition probability on lv18


```R
#path length matrix
pathLength <- matrix(0,5,5)
rownames(pathLength) <- c(3,2,4,5,1)
colnames(pathLength) <- c(3,2,4,5,1)
for(i in c(3,2,4,5,1)){
  for(j in c(3,2,4,5,1)){
    MarkovHCPath <- findTransitionPath(MarkovObject = MarkovHC_dt,
                                       level = 18,
                                       basinA = i,
                                       basinB = j)
    if(length(MarkovHCPath)>0){
      pathLength[which(rownames(pathLength)==i),which(colnames(pathLength)==j)] <- MarkovHCPath[[3]]
    }else{
      pathLength[which(rownames(pathLength)==i),which(colnames(pathLength)==j)] <- Inf
    }
  }
}
```


```R
pathLength
```


<table>
<caption>A matrix: 5 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>3</th><th scope=col>2</th><th scope=col>4</th><th scope=col>5</th><th scope=col>1</th></tr>
</thead>
<tbody>
	<tr><th scope=row>3</th><td>0.000000</td><td>2.906242</td><td>3.801976</td><td>5.922213</td><td>10.771162</td></tr>
	<tr><th scope=row>2</th><td>1.883931</td><td>0.000000</td><td>3.981589</td><td>3.015971</td><td>10.534989</td></tr>
	<tr><th scope=row>4</th><td>2.840507</td><td>4.927288</td><td>0.000000</td><td>5.361802</td><td> 8.549278</td></tr>
	<tr><th scope=row>5</th><td>4.868345</td><td>2.855061</td><td>5.535095</td><td>0.000000</td><td> 7.519018</td></tr>
	<tr><th scope=row>1</th><td>2.433467</td><td>4.800132</td><td>2.924561</td><td>1.945071</td><td> 0.000000</td></tr>
</tbody>
</table>




```R
P_matrix <- pathLength

P_matrix <- exp(-P_matrix*2)

P_matrix
```


<table>
<caption>A matrix: 5 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>3</th><th scope=col>2</th><th scope=col>4</th><th scope=col>5</th><th scope=col>1</th></tr>
</thead>
<tbody>
	<tr><th scope=row>3</th><td>1.0000000000</td><td>0.0029899941</td><td>4.984771e-04</td><td>7.178460e-06</td><td>4.408469e-10</td></tr>
	<tr><th scope=row>2</th><td>0.0231013990</td><td>1.0000000000</td><td>3.480454e-04</td><td>2.400828e-03</td><td>7.070086e-10</td></tr>
	<tr><th scope=row>4</th><td>0.0034100958</td><td>0.0000525064</td><td>1.000000e+00</td><td>2.201900e-05</td><td>3.751386e-08</td></tr>
	<tr><th scope=row>5</th><td>0.0000590758</td><td>0.0033122716</td><td>1.556962e-05</td><td>1.000000e+00</td><td>2.944854e-07</td></tr>
	<tr><th scope=row>1</th><td>0.0076969228</td><td>0.0000677109</td><td>2.882431e-03</td><td>2.044244e-02</td><td>1.000000e+00</td></tr>
</tbody>
</table>




```R
diag(P_matrix) <- 0
```


```R
P_matrix <- P_matrix/rowSums(P_matrix)
```


```R
P_matrix
```


<table>
<caption>A matrix: 5 × 5 of type dbl</caption>
<thead>
	<tr><th></th><th scope=col>3</th><th scope=col>2</th><th scope=col>4</th><th scope=col>5</th><th scope=col>1</th></tr>
</thead>
<tbody>
	<tr><th scope=row>3</th><td>0.00000000</td><td>0.855347090</td><td>0.142599243</td><td>0.002053541</td><td>1.261130e-07</td></tr>
	<tr><th scope=row>2</th><td>0.89366172</td><td>0.000000000</td><td>0.013463896</td><td>0.092874359</td><td>2.735014e-08</td></tr>
	<tr><th scope=row>4</th><td>0.97860252</td><td>0.015067876</td><td>0.000000000</td><td>0.006318840</td><td>1.076543e-05</td></tr>
	<tr><th scope=row>5</th><td>0.01744084</td><td>0.977875634</td><td>0.004596589</td><td>0.000000000</td><td>8.694035e-05</td></tr>
	<tr><th scope=row>1</th><td>0.24757305</td><td>0.002177935</td><td>0.092713977</td><td>0.657535040</td><td>0.000000e+00</td></tr>
</tbody>
</table>




```R
pheatmap::pheatmap(as.matrix(P_matrix), cluster_rows = F, cluster_cols =F,
                   scale = "none" ,
                   legend_breaks= ceiling(seq(min(P_matrix),
                                              1,0.01)),
                   color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(P_matrix),1,0.01))),
                   breaks= seq(min(P_matrix),
                               1,
                               by=0.01),
                   show_colnames = T, show_rownames = T,
                   #annotation_col  = annotation_col_C_andcluster,
                   #annotation_colors = ann_colors_C,
                   width=10, 
                   heigheight = 10,
                   fontsize =10,
                   cellwidth = 30,
                   cellheight = 30,
                   display_numbers=TRUE,
                   number_color = 'black',
                   number_format = "%.4f",
                   filename = './PmatrixLv18.pdf'
)
```


```R
options(repr.plot.width=5, repr.plot.height=5)
scales::show_col(c("#08519c",
                   "#238b45",
                   "#cb181d",
                   "#9ecae1",
                   "#fd8d3c"))
```


![png](output_48_0.png)



```R
layout$lv18 <- MarkovHCLabels$lv18
```


```R
ggplot(data=layout, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=lv18), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
  scale_fill_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5'))
```


![png](output_50_0.png)



```R
pdf(file = './lv18.pdf', width = 3.5, height = 3.5)
ggplot(data=layout, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=lv18), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
  scale_fill_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5'))
dev.off()
```


<strong>png:</strong> 2


# heatmap


```R
Idents(dt_sub) <- MarkovHCLabels$lv18
```


```R
markers <- FindAllMarkers(dt_sub,
                        min.pct = 0.25,
                        logfc.threshold = 0.25,
                        only.pos=TRUE)
```

    Calculating cluster 4
    
    Calculating cluster 2
    
    Calculating cluster 5
    
    Calculating cluster 3
    
    Calculating cluster 1
    



```R
markerstop50 <- markers %>% group_by(cluster) %>% top_n(n = 50, wt = avg_logFC)%>%as.data.frame()
```


```R
head(markerstop50)
```


<table>
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>1</th><td>2.514655e-137</td><td>0.9800802</td><td>0.997</td><td>1</td><td>5.761074e-133</td><td>4</td><td>RPL34  </td></tr>
	<tr><th scope=row>2</th><td>2.591188e-136</td><td>1.1945679</td><td>0.997</td><td>1</td><td>5.936411e-132</td><td>4</td><td>RPL7   </td></tr>
	<tr><th scope=row>3</th><td>4.364072e-123</td><td>0.8147982</td><td>1.000</td><td>1</td><td>9.998089e-119</td><td>4</td><td>MT-CO2 </td></tr>
	<tr><th scope=row>4</th><td>5.428021e-109</td><td>0.8003757</td><td>0.997</td><td>1</td><td>1.243560e-104</td><td>4</td><td>RPS3A  </td></tr>
	<tr><th scope=row>5</th><td>8.623314e-109</td><td>0.6452350</td><td>1.000</td><td>1</td><td>1.975601e-104</td><td>4</td><td>MT-ATP6</td></tr>
	<tr><th scope=row>6</th><td>3.586057e-103</td><td>0.7746922</td><td>1.000</td><td>1</td><td> 8.215657e-99</td><td>4</td><td>RPS27A </td></tr>
</tbody>
</table>




```R
Idents(dt_sub) <- factor(Idents(dt_sub), levels = c('3','2','4','5','1'))
```


```R
dt_sub@meta.data$cellTypes <- mapvalues(as.character(rownames(dt_sub@meta.data)), from=as.character(rownames(dt_update@meta.data)), to=as.character(Idents(dt_update)), warn_missing = FALSE)
```


```R
#heatmap
ordergenes_expression_matrix <- GetAssayData(object = dt_sub, slot = "scale.data")%>%as.data.frame()
markerstop50 <- markerstop50[order(factor(markerstop50$cluster, levels =c('3','2','4','5','1')), decreasing = FALSE),]
ordergenes_expression_matrix <- subset(ordergenes_expression_matrix, rownames(ordergenes_expression_matrix)%in%markerstop50$gene)
ordergenes_expression_matrix <- ordergenes_expression_matrix[order(factor(rownames(ordergenes_expression_matrix), levels = unique(markerstop50$gene))),]

annotation_col_C_andcluster = data.frame(Basins=factor(dt_sub@meta.data$basins),
                                         cellTypes=factor(dt_sub@meta.data$cellTypes))
rownames(annotation_col_C_andcluster) = colnames(ordergenes_expression_matrix)
ann_colors_C = list(
  Basins =  c("1"=alpha("#fd8d3c",1),
              "2"=alpha("#238b45",1),
              "3"=alpha("#08519c",1),
              "4"=alpha("#9ecae1",1),
              "5"=alpha("#cb181d",1)),
  cellTypes = c("Cancer cell"=alpha("#cb181d",1),
                "MSC"=alpha("#08519c",1)))
```


```R
write.table(markerstop50, file = './markerGenes.txt', quote = FALSE, sep=' ', row.names = F)
```

# write marker genes


```R
head(markerstop50)
```


<table>
<caption>A data.frame: 6 × 8</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th><th scope=col>basins</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>151</th><td>4.062769e-41</td><td>0.6867383</td><td>0.982</td><td>0.909</td><td>9.307803e-37</td><td>3</td><td>CLDN4   </td><td>3</td></tr>
	<tr><th scope=row>152</th><td>2.867994e-40</td><td>0.8157101</td><td>0.763</td><td>0.545</td><td>6.570573e-36</td><td>3</td><td>EDN1    </td><td>3</td></tr>
	<tr><th scope=row>153</th><td>2.158714e-37</td><td>0.7612563</td><td>0.966</td><td>0.887</td><td>4.945613e-33</td><td>3</td><td>UBC     </td><td>3</td></tr>
	<tr><th scope=row>154</th><td>3.427924e-33</td><td>0.5055630</td><td>0.948</td><td>0.821</td><td>7.853375e-29</td><td>3</td><td>DDX5    </td><td>3</td></tr>
	<tr><th scope=row>155</th><td>2.243492e-31</td><td>0.4927032</td><td>0.895</td><td>0.806</td><td>5.139840e-27</td><td>3</td><td>CIRBP   </td><td>3</td></tr>
	<tr><th scope=row>156</th><td>4.292957e-31</td><td>0.5829897</td><td>0.708</td><td>0.475</td><td>9.835163e-27</td><td>3</td><td>PDZK1IP1</td><td>3</td></tr>
</tbody>
</table>




```R
markerstop50$basins <- as.character(markerstop50$cluster)
```


```R
markerstop50$basins <- mapvalues(markerstop50$basins, 
                                                from=c('3','2','4','5','1'), 
                                                to=c('basin1','basin2','basin3','basin4','basin5'))
```


```R
write.table(markerstop50, file = './markerGenes.txt', quote = FALSE, sep=',', row.names = F)
```


```R
ordered_genes_expression_matrix <- t(ordergenes_expression_matrix)%>%as.data.frame()
ordered_genes_expression_matrix$Basins <- factor(dt_sub@meta.data$basins, levels =c('3','2','4','5','1'))
ordered_genes_expression_matrix$cellTypes <- factor(dt_sub@meta.data$cellTypes, levels = rev(c( "MSC","Cancer cell")))
```


```R
dim(ordered_genes_expression_matrix)
```


<style>
.list-inline {list-style: none; margin:0; padding: 0}
.list-inline>li {display: inline-block}
.list-inline>li:not(:last-child)::after {content: "\00b7"; padding: 0 .5ex}
</style>
<ol class=list-inline><li>1526</li><li>239</li></ol>




```R
ordered_genes_expression_matrix <- doBy::orderBy(~ Basins + cellTypes, ordered_genes_expression_matrix)
```


```R
ordered_genes_expression_matrix <- ordered_genes_expression_matrix[,-which(colnames(ordered_genes_expression_matrix)%in%c('Basins','cellTypes'))]
ordered_genes_expression_matrix <- t(ordered_genes_expression_matrix)
ordered_genes_expression_matrix_copy <- ordered_genes_expression_matrix
```


```R
ordered_genes_expression_matrix[ordered_genes_expression_matrix>4] <- 4
ordered_genes_expression_matrix[ordered_genes_expression_matrix< (-4)] <- (-4)
pheatmap(as.matrix(ordered_genes_expression_matrix), cluster_rows = F, cluster_cols =F,
         scale = "none" ,
         legend_breaks= ceiling(seq(min(ordered_genes_expression_matrix),
                                    max(ordered_genes_expression_matrix),0.01)),
         color = colorRampPalette(colors = c("#377eb8","#deebf7","#e41a1c"))(length(seq(min(ordered_genes_expression_matrix),max(ordered_genes_expression_matrix),0.01))),
         breaks= seq(min(ordered_genes_expression_matrix),
                     max(ordered_genes_expression_matrix),
                     by=0.01),
         show_colnames = F, show_rownames = T,
         annotation_col  = annotation_col_C_andcluster,
         annotation_colors = ann_colors_C,width=5, heigheight = 180,
         fontsize =2,
         filename = './markerGeneHeatmap.pdf'
         )
```

# GO


```R
head(markerstop50,n=3)
```


<table>
<caption>A data.frame: 3 × 7</caption>
<thead>
	<tr><th></th><th scope=col>p_val</th><th scope=col>avg_logFC</th><th scope=col>pct.1</th><th scope=col>pct.2</th><th scope=col>p_val_adj</th><th scope=col>cluster</th><th scope=col>gene</th></tr>
	<tr><th></th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;chr&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>151</th><td>4.062769e-41</td><td>0.6867383</td><td>0.982</td><td>0.909</td><td>9.307803e-37</td><td>3</td><td>CLDN4</td></tr>
	<tr><th scope=row>152</th><td>2.867994e-40</td><td>0.8157101</td><td>0.763</td><td>0.545</td><td>6.570573e-36</td><td>3</td><td>EDN1 </td></tr>
	<tr><th scope=row>153</th><td>2.158714e-37</td><td>0.7612563</td><td>0.966</td><td>0.887</td><td>4.945613e-33</td><td>3</td><td>UBC  </td></tr>
</tbody>
</table>




```R
markerstop50$cluster <- as.character(markerstop50$cluster)
basins <- unique(markerstop50[,6])
for(i in 1:length(basins)){
  upregulatedGenes <- (subset(markerstop50, markerstop50[,6]==basins[i])%>%as.data.frame())[,7]
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
  write.table(GO_upregulatedGenes.result, file = paste('./GO/',as.character(basins[i]),'.txt', sep=''))
}
```

# lv21


```R
layout$lv21 <- MarkovHCLabels$lv21
```


```R
pdf(file = './lv21.pdf', width = 3.5, height = 3.5)
ggplot(data=layout, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=lv21), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
  scale_fill_manual(
    values = c( "2"=alpha("#cb181d",0.7),
                "1"=alpha("#08519c",0.7)),
    breaks = c( "1",
                 "2"))
dev.off()
```


<strong>png:</strong> 2


# feature plot


```R
# MSC, the intestinal stem cell markers： "OLFM4", "EPHB2", "SOX9"
# cancer : CEACAM5, CEACAM6, BAX, CCND2
# GC-related early diagnosis markers: FABP1, CEACAM5, and CDH17
# EGC panel: DPCR1, MUC5AC, KLK10, SLC11A2, SULT2B1, KLK7, ECM1, LMTK3
options(repr.plot.width=9, repr.plot.height=8)
pdf('./markerGenes.pdf', height = 9, width = 10)
FeaturePlot(dt_sub, features = c("OLFM4", "EPHB2", "SOX9",
                                  "DPCR1","MUC5AC", "KLK10",  "ECM1", "SULT2B1", "LMTK3"))
dev.off()
```


<strong>png:</strong> 2


# transition paths


```R
stepWisepath1 = stepWisepath(MarkovObject=MarkovHC_dt,
                        MarkovLevel=18,
                        stepBasin=c(3,2,5))
```


```R
stepWisepath2 = stepWisepath(MarkovObject=MarkovHC_dt,
                        MarkovLevel=18,
                        stepBasin=c(3,4,5))
```


```R
basins <- MarkovHCLabels$lv18
```


```R
basins[stepWisepath1[[2]][[2]]] <- names(stepWisepath1[[2]])[2]

basins[stepWisepath1[[2]][[1]]] <- names(stepWisepath1[[2]])[1]

basins[stepWisepath2[[2]][[2]]] <- names(stepWisepath2[[2]])[2]

basins[stepWisepath2[[2]][[1]]] <- names(stepWisepath2[[2]])[1]
```


```R
options(repr.plot.width=5, repr.plot.height=5)
dt_sub@meta.data$basins <- basins
DimPlot(dt_sub, reduction = "umap", group.by = 'basins',label=T,pt.size=2, label.size=5)
```


![png](output_84_0.png)



```R
layout$detail <- dt_sub@meta.data$basins
```


```R
bp <- subset(layout, layout$detail %in% c('1', '2', '3', '4','5'))
cp <- subset(layout, layout$detail %in% c(names(stepWisepath1[[2]])[1],names(stepWisepath1[[2]])[2],names(stepWisepath2[[2]])[1],names(stepWisepath2[[2]])[2]))
```


```R
ggplot(data=bp, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=detail), color=alpha("#525252",0))+
  geom_point(data=cp,size=1, shape=21, aes(x=UMAP_1, y=UMAP_2, fill=detail), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
   scale_fill_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))
```


![png](output_87_0.png)



```R
pdf(file = './detailnolegend.pdf', width = 3.5, height = 3.5)
ggplot(data=bp, mapping =  aes(x=UMAP_1, y=UMAP_2)) +
  geom_point(size=1, shape=21, aes(fill=detail), color=alpha("#525252",0))+
  geom_point(data=cp,size=1, shape=21, aes(x=UMAP_1, y=UMAP_2, fill=detail), color=alpha("#525252",0))+
  xlim(min(layout$UMAP_1)-1,max(layout$UMAP_1)+1)+
  ylim(min(layout$UMAP_2)-1,max(layout$UMAP_2)+1)+guides(fill=FALSE)+
  mytheme+
  xlab("UMAP_1") + ylab("UMAP_2")+
   scale_fill_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))
dev.off()
```


<strong>png:</strong> 2


# path1 heatmap


```R
path1 <- c(rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '3')], 
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == 'cp32')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '2')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == 'cp25')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '5')])
```


```R
#get counts
dt_counts <- GetAssayData(object = dt_sub, slot = "counts")%>%as.matrix()%>%t()
dt_counts <- subset(dt_counts, rownames(dt_counts)%in%path1)
dt_counts <- dt_counts[order(factor(rownames(dt_counts), levels = unique(path1)), decreasing = FALSE),]
dt_counts <- as.data.frame(dt_counts)
pseudotime <- as.data.frame(1:nrow(dt_counts))
```


```R
#monocle object
gene_metadata <- as.data.frame(colnames(dt_counts))
rownames(gene_metadata) <- gene_metadata[,1]
gene_metadata $ gene_short_name <- gene_metadata[,1]
colnames(gene_metadata) <- c('gene_short_name','ensembleID')
rownames(pseudotime) <- rownames(dt_counts)

path_object <- newCellDataSet(as.matrix(t(dt_counts)),
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

path_diff_test_res <- differentialGeneTest(path_object[dt_sub@assays$RNA@var.features,],
                                           fullModelFormulaStr = "~sm.ns(Pseudotime)")
path_sig_gene_names <- row.names(subset(path_diff_test_res, qval < 0.1))
```

    Removing 20 outliers
    



```R
path_pheatmap <- plot_pseudotime_heatmap(path_object[path_sig_gene_names,],
                                         num_clusters = 6,
                                         cores = 1,
                                         show_rownames = F,
                                         return_heatmap = TRUE)
```

    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>



![png](output_93_1.png)



```R
geneclusters <- cutree(path_pheatmap$tree_row, k=6)
orderedGenes <- c()
for(i in c(2,3,5,6,4,1)){
  orderedGenes <- c(orderedGenes, names(geneclusters[which(geneclusters==i)]))
}
```


```R
path_diff_test_res_write <- path_diff_test_res[rownames(path_diff_test_res)%in%orderedGenes,]
```


```R
path_diff_test_res_write$gene_short_name <- factor(path_diff_test_res_write$gene_short_name, levels = orderedGenes)
path_diff_test_res_write <- path_diff_test_res_write[order(path_diff_test_res_write$gene_short_name, decreasing = FALSE),]
```


```R
head(path_diff_test_res_write)
```


<table>
<caption>A data.frame: 6 × 7</caption>
<thead>
	<tr><th></th><th scope=col>status</th><th scope=col>family</th><th scope=col>pval</th><th scope=col>qval</th><th scope=col>gene_short_name</th><th scope=col>ensembleID</th><th scope=col>num_cells_expressed</th></tr>
	<tr><th></th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;dbl&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;fct&gt;</th><th scope=col>&lt;dbl&gt;</th></tr>
</thead>
<tbody>
	<tr><th scope=row>MMP7</th><td>OK</td><td>negbinomial.size</td><td> 1.259013e-47</td><td> 3.342512e-46</td><td>MMP7 </td><td>MMP7 </td><td> 99</td></tr>
	<tr><th scope=row>PGC</th><td>OK</td><td>negbinomial.size</td><td>3.046339e-167</td><td>3.515006e-165</td><td>PGC  </td><td>PGC  </td><td> 99</td></tr>
	<tr><th scope=row>DPCR1</th><td>OK</td><td>negbinomial.size</td><td> 7.192836e-03</td><td> 1.542424e-02</td><td>DPCR1</td><td>DPCR1</td><td> 13</td></tr>
	<tr><th scope=row>MSMB</th><td>OK</td><td>negbinomial.size</td><td> 8.508578e-06</td><td> 2.907259e-05</td><td>MSMB </td><td>MSMB </td><td> 62</td></tr>
	<tr><th scope=row>LCN2</th><td>OK</td><td>negbinomial.size</td><td> 1.770022e-09</td><td> 9.108177e-09</td><td>LCN2 </td><td>LCN2 </td><td>420</td></tr>
	<tr><th scope=row>TFF3</th><td>OK</td><td>negbinomial.size</td><td> 9.365198e-65</td><td> 3.796702e-63</td><td>TFF3 </td><td>TFF3 </td><td>850</td></tr>
</tbody>
</table>




```R
write.table(path_diff_test_res_write, file = './path_diff_test_res_writepath1.txt', quote = F, sep = ' ')
```


```R
Path1.pdf <- plot_pseudotime_heatmap(path_object[orderedGenes,],
                        num_clusters = 1,
                        cores = 1,
                        show_rownames = T,
                        cluster_rows=FALSE,
                         return_heatmap = TRUE)
```

    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>



![png](output_99_1.png)



```R
pdf("./Path1.pdf", height = 15, width = 6)
Path1.pdf
dev.off()
```


<strong>png:</strong> 2



```R
dt_scaledData <- GetAssayData(object = dt_sub, slot = "scale.data")%>%as.matrix()%>%t()
dt_scaledData <- subset(dt_scaledData, rownames(dt_scaledData)%in%path1)
dt_scaledData <- dt_scaledData[order(factor(rownames(dt_scaledData), levels = unique(path1)), decreasing = FALSE),]
dt_scaledData <- as.data.frame(dt_scaledData)
```


```R
dt_scaledData$pseudotime <- 1:nrow(dt_scaledData)
```


```R
dt_scaledData$basins <- dt_sub@meta.data[path1,]$basins
```


```R
ggplot(dt_scaledData,
       aes(x=pseudotime,y=OLFM4,color=basins)) + geom_point(shape=19, size=2.5) +
  ylab("scaled data") + xlab("pseudo time") + ggtitle('OLFM4')+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+
   scale_colour_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))+mytheme
```

    `geom_smooth()` using formula 'y ~ x'
    



![png](output_104_1.png)



```R
for(i in c('OLFM4', 'CEACAM6')){

pdf(paste('./path1_',i,'.pdf', sep=''),width = 3.5, height = 3.5)
print(ggplot(dt_scaledData,
       aes(x=pseudotime,y=get(i),color=basins)) + geom_point(shape=19, size=1) +
  ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+
   scale_colour_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))+mytheme)
dev.off()
    
}
```

    `geom_smooth()` using formula 'y ~ x'
    
    `geom_smooth()` using formula 'y ~ x'
    


# differentially expressed genes for GO


```R
geneclusters <- cutree(path_pheatmap$tree_row, k=6)
increasedGenes <- c()
for(i in c(6,4,1)){
  increasedGenes <- c(increasedGenes, names(geneclusters[which(geneclusters==i)]))
}
increasedGenes_path1 <- increasedGenes

decreasedGenes <- c()
for(i in c(2,3,5)){
  decreasedGenes <- c(decreasedGenes, names(geneclusters[which(geneclusters==i)]))
}
decreasedGenes_path1 <- decreasedGenes
```


```R
GO_increasedGenes <- enrichGO(gene = increasedGenes,
                              keyType = "SYMBOL",
                              OrgDb = 'org.Hs.eg.db',
                              ont = "BP",
                              pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.2,
                              minGSSize = 3,
                              maxGSSize = 500,
                              readable = FALSE)
GO_increasedGenes.result <- as.data.frame(GO_increasedGenes@result)
write.table(GO_increasedGenes.result, file = './GO_increasedGenes.path1.txt', sep=',')
```


```R
GO_decreasedGenes <- enrichGO(gene = decreasedGenes,
                              keyType = "SYMBOL",
                              OrgDb = 'org.Hs.eg.db',
                              ont = "BP",
                              pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.2,
                              minGSSize = 3,
                              maxGSSize = 500,
                              readable = FALSE)
GO_decreasedGenes.result <- as.data.frame(GO_decreasedGenes@result)
write.table(GO_decreasedGenes.result, file = './GO_decreasedGenes.path1.txt', sep=',')
```

# path2 heatmap


```R
path2 <- c(rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '3')], 
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == 'cp34')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '4')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == 'cp45')],
           rownames(dt_sub@meta.data)[which(dt_sub@meta.data$basins == '5')])
```


```R
#get counts
dt_counts <- GetAssayData(object = dt_sub, slot = "counts")%>%as.matrix()%>%t()
dt_counts <- subset(dt_counts, rownames(dt_counts)%in%path2)
dt_counts <- dt_counts[order(factor(rownames(dt_counts), levels = unique(path2)), decreasing = FALSE),]
dt_counts <- as.data.frame(dt_counts)
pseudotime <- as.data.frame(1:nrow(dt_counts))
```


```R
#monocle object
gene_metadata <- as.data.frame(colnames(dt_counts))
rownames(gene_metadata) <- gene_metadata[,1]
gene_metadata $ gene_short_name <- gene_metadata[,1]
colnames(gene_metadata) <- c('gene_short_name','ensembleID')
rownames(pseudotime) <- rownames(dt_counts)

path_object <- newCellDataSet(as.matrix(t(dt_counts)),
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

path_diff_test_res <- differentialGeneTest(path_object[dt_sub@assays$RNA@var.features,],
                                           fullModelFormulaStr = "~sm.ns(Pseudotime)")
path_sig_gene_names <- row.names(subset(path_diff_test_res, qval < 0.1))
```

    Removing 19 outliers
    



```R
path_pheatmap <- plot_pseudotime_heatmap(path_object[path_sig_gene_names,],
                                         num_clusters = 5,
                                         cores = 1,
                                         show_rownames = F,
                                         return_heatmap = TRUE)
```

    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>



![png](output_114_1.png)



```R
geneclusters <- cutree(path_pheatmap$tree_row, k=5)
orderedGenes <- c()
for(i in c(3,4,5,2,1)){
  orderedGenes <- c(orderedGenes, names(geneclusters[which(geneclusters==i)]))
}
```


```R
path_diff_test_res_write <- path_diff_test_res[rownames(path_diff_test_res)%in%orderedGenes,]
path_diff_test_res_write$gene_short_name <- factor(path_diff_test_res_write$gene_short_name, levels = orderedGenes)
path_diff_test_res_write <- path_diff_test_res_write[order(path_diff_test_res_write$gene_short_name, decreasing = FALSE),]
```


```R
write.table(path_diff_test_res_write, file = './path_diff_test_res_writepath2.txt', quote = F, sep = ' ')
```


```R
Path2.pdf <- plot_pseudotime_heatmap(path_object[orderedGenes,],
                                            num_clusters = 1,
                                            cores = 1,
                                            show_rownames = T,
                                            cluster_rows=FALSE,
                                            return_heatmap = TRUE)
```

    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>
    <simpleError in checkwz(wz, M = M, trace = trace, wzepsilon = control$wzepsilon): NAs found in the working weights variable 'wz'>



![png](output_118_1.png)



```R
pdf("./Path2.pdf", height = 15, width = 6)
Path2.pdf
dev.off()
```


<strong>png:</strong> 2



```R
dt_scaledData <- GetAssayData(object = dt_sub, slot = "scale.data")%>%as.matrix()%>%t()
dt_scaledData <- subset(dt_scaledData, rownames(dt_scaledData)%in%path2)
dt_scaledData <- dt_scaledData[order(factor(rownames(dt_scaledData), levels = unique(path2)), decreasing = FALSE),]
dt_scaledData <- as.data.frame(dt_scaledData)
```


```R
dt_scaledData$pseudotime <- 1:nrow(dt_scaledData)
```


```R
dt_scaledData$basins <- dt_sub@meta.data[path2,]$basins
```


```R
ggplot(dt_scaledData,
       aes(x=pseudotime,y=SOX4,color=basins)) + geom_point(shape=19, size=2.5) +
  ylab("scaled data") + xlab("pseudo time") + ggtitle('SOX4')+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+
   scale_colour_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))+mytheme
```

    `geom_smooth()` using formula 'y ~ x'
    



![png](output_123_1.png)



```R
for(i in c('SOX4', 'NEAT1')){

pdf(paste('./path2_',i,'.pdf', sep=''),width = 3.5, height = 3.5)
print(ggplot(dt_scaledData,
       aes(x=pseudotime,y=get(i),color=basins)) + geom_point(shape=19, size=1) +
  ylab("scaled data") + xlab("pseudo time") + ggtitle(i)+
  geom_smooth(method = loess,fill='#969696',colour='#737373')+
   scale_colour_manual(
    values = c( "1"=alpha("#fd8d3c",0.7),
                "2"=alpha("#238b45",0.7),
                "3"=alpha("#08519c",0.7),
                "4"=alpha("#9ecae1",0.7),
                "5"=alpha("#cb181d",0.7),
               'cp32'=alpha('#54278f',0.7),
               'cp25'=alpha('#54278f',0.7),
               'cp34'=alpha('#54278f',0.7),
               'cp45'=alpha('#54278f',0.7)),
    breaks = c( "1",
                "2",
                "3",
                '4',
                '5',
                'cp32',
                'cp25',
                'cp34',
                'cp45'))+mytheme)
dev.off()
    
}
```

    `geom_smooth()` using formula 'y ~ x'
    
    `geom_smooth()` using formula 'y ~ x'
    


# differentially expressed genes for GO


```R
geneclusters <- cutree(path_pheatmap$tree_row, k=5)
increasedGenes <- c()
for(i in c(1,2)){
  increasedGenes <- c(increasedGenes, names(geneclusters[which(geneclusters==i)]))
}
increasedGenes_path253 <- increasedGenes

decreasedGenes <- c()
for(i in c(3,4,5)){
  decreasedGenes <- c(decreasedGenes, names(geneclusters[which(geneclusters==i)]))
}
decreasedGenes_path253 <- decreasedGenes
```


```R
GO_increasedGenes <- enrichGO(gene = increasedGenes,
                              keyType = "SYMBOL",
                              OrgDb = 'org.Hs.eg.db',
                              ont = "BP",
                              pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.2,
                              minGSSize = 3,
                              maxGSSize = 500,
                              readable = FALSE)
GO_increasedGenes.result <- as.data.frame(GO_increasedGenes@result)
write.table(GO_increasedGenes.result, file = './GO_increasedGenes.path2.txt', sep=',')

GO_decreasedGenes <- enrichGO(gene = decreasedGenes,
                              keyType = "SYMBOL",
                              OrgDb = 'org.Hs.eg.db',
                              ont = "BP",
                              pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.2,
                              minGSSize = 3,
                              maxGSSize = 500,
                              readable = FALSE)
GO_decreasedGenes.result <- as.data.frame(GO_decreasedGenes@result)
write.table(GO_decreasedGenes.result, file = './GO_decreasedGenes.path2.txt', sep=',')
```

# step wise DEG


```R
Idents(dt_sub) <- MarkovHCLabels$lv18
```


```R
paiwise <- matrix(c(3,2,3,4,2,5,4,5), ncol = 2)
```


```R
for(i in 1:nrow(paiwise)){
    Marker_temp <- FindMarkers(dt_sub, 
                               only.pos = FALSE,
                               test.use = 'wilcox',
                               min.pct = 0.1, 
                               logfc.threshold = 0.1,
                               ident.1 = paiwise[i,1], 
                               ident.2 = paiwise[i,2])
    
    Marker_temp_up <- subset(Marker_temp, Marker_temp$avg_logFC>0)
    GO_upGenes <- enrichGO(gene = rownames(Marker_temp_up),
                              keyType = "SYMBOL",
                              OrgDb = 'org.Hs.eg.db',
                              ont = "ALL",
                              pAdjustMethod = "fdr",
                              pvalueCutoff = 0.05,
                              qvalueCutoff  = 0.2,
                              minGSSize = 3,
                              maxGSSize = 500,
                              readable = FALSE)
   GO_upGenes.result <- as.data.frame(GO_upGenes@result)
   write.table(GO_upGenes.result, file =paste('./stepwiseDEG/GO_upGenes.',as.character(paiwise[i,1]),'vs',as.character(paiwise[i,2]),'.txt', sep=''))

   Marker_temp_down <- subset(Marker_temp, Marker_temp$avg_logFC<0)
   GO_downGenes <- enrichGO(gene = rownames(Marker_temp_down),
                       keyType = "SYMBOL",
                       OrgDb = 'org.Hs.eg.db',
                       ont = "ALL",
                       pAdjustMethod = "fdr",
                       pvalueCutoff = 0.05,
                       qvalueCutoff  = 0.2,
                       minGSSize = 3,
                       maxGSSize = 500,
                       readable = FALSE)
  GO_downGenes.result <- as.data.frame(GO_downGenes@result)
  write.table(GO_downGenes.result, file =paste('./stepwiseDEG/GO_downGenes.',as.character(paiwise[i,1]),'vs',as.character(paiwise[i,2]),'.txt', sep=''))
}
```


```R
save.image('./gastric_cancer.RData')
```
