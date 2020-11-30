```R
library(dplyr)
library(car)
library(dplyr)
library(MASS)
library(MarkovHC)
library(ggplot2)
library(ggalluvial)
library(ggraph)
library(viridis)
library(stringr)
library(ggforce)
library(doBy)
library(alluvial)
library(rgl)
library(plot3D)
```

    
    Attaching package: ‘dplyr’
    
    
    The following objects are masked from ‘package:stats’:
    
        filter, lag
    
    
    The following objects are masked from ‘package:base’:
    
        intersect, setdiff, setequal, union
    
    
    Loading required package: carData
    
    
    Attaching package: ‘car’
    
    
    The following object is masked from ‘package:dplyr’:
    
        recode
    
    
    
    Attaching package: ‘MASS’
    
    
    The following object is masked from ‘package:dplyr’:
    
        select
    
    
    Loading required package: parallel
    
    Loading required package: doParallel
    
    Loading required package: foreach
    
    Loading required package: iterators
    
    Loading required package: dbscan
    
    Loading required package: igraph
    
    
    Attaching package: ‘igraph’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        as_data_frame, groups, union
    
    
    The following objects are masked from ‘package:stats’:
    
        decompose, spectrum
    
    
    The following object is masked from ‘package:base’:
    
        union
    
    
    Loading required package: Matrix
    
    Loading required package: Rcpp
    
    Loading required package: plyr
    
    ------------------------------------------------------------------------------
    
    You have loaded plyr after dplyr - this is likely to cause problems.
    If you need functions from both plyr and dplyr, please load plyr first, then dplyr:
    library(plyr); library(dplyr)
    
    ------------------------------------------------------------------------------
    
    
    Attaching package: ‘plyr’
    
    
    The following objects are masked from ‘package:dplyr’:
    
        arrange, count, desc, failwith, id, mutate, rename, summarise,
        summarize
    
    
    Loading required package: doBy
    
    Loading required package: ggraph
    
    Loading required package: ggplot2
    
    Loading required package: viridisLite
    



    Error in library(alluvial): there is no package called ‘alluvial’
    Traceback:


    1. library(alluvial)



```R
data("example_3D1",package="MarkovHC")
```


```R
example5 <- example_3D1
```

# run MarkovHC


```R
MarkovHC_5 <- MarkovHC(origin_matrix=t(example5),
                       transformtype="none",
                       KNN=40,
                       basecluster="kmeans",
                       dobasecluster=FALSE,
                       baseclusternum=100,
                       emphasizedistance=1,
                       weightDist=2,
                       weightDens=0.5,
                       cutpoint=0.01,
                       showprocess=FALSE,
                       bn=2,
                       minBasinSize=0.2,
                       noiseBasinSize=5)
```


```R
labels <-  fetchLabels(MarkovObject=MarkovHC_5,
                       MarkovLevels=1:length(MarkovHC_5$hierarchicalStructure))

basin <- labels[,9]

for(i in 1:length(basin)){
    basin[i] <- str_split( basin[i], '\\+')[[1]][1]
}
```


```R
basin <- paste('basin',basin,sep = '')
level <- 9
for (i in 1:length(MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
```


```R
#basins
colors <- c("#e41a1c","#377eb8","#4daf4a")
colors <- colors[factor(basin)]

scatter3d(x=example5[,1],
          y=example5[,2], 
          z=example5[,3],
          point.col=colors,
          pch = 21,
          sphere.size=1.5,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
view3d(userMatrix = um)
rgl.postscript("./TransitionGroup.3dsimple.basins.pdf", "pdf", drawText = FALSE)

#find the transition points
label <- findTransitionPoints(MarkovObject = MarkovHC_5,
                              level = level,
                              basinA = 1,
                              basinB = 2)
basin[which(label==1)] <- 'transition points\nbetween 1 and 2'

colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3")
colors <- colors[factor(basin)]

scatter3d(x=example5[,1],
          y=example5[,2], 
          z=example5[,3],
          point.col=colors,
          pch = 21,
          sphere.size=1.5,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
um <- par3d()$userMatrix
view3d(userMatrix = um)
rgl.snapshot("./TransitionGroup.3dsimple.1.png", fmt = "png", top = TRUE )
rgl.postscript("./TransitionGroup.3dsimple.1.pdf", "pdf", drawText = FALSE)

#find the transition paths
level <- 9
label <- findTransitionPath(MarkovObject = MarkovHC_5,
                            level = level,
                            basinA = 1,
                            basinB = 2)
basin[which(label==1)] <- 'transition path\nfrom 1 to 2'

for (i in 1:length(MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example5[,1],
          y=example5[,2], 
          z=example5[,3],
          point.col=colors,
          pch = 21,
          sphere.size=1.5,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
view3d(userMatrix = um)
rgl.postscript("./TransitionGroup.3dsimple.2.pdf", "pdf", drawText = FALSE)

level <- 9
label <- findTransitionPath(MarkovObject = MarkovHC_5,
                            level = level,
                            basinA = 2,
                            basinB = 1)
basin[which(label==1)] <- 'transition path\nfrom 2 to 1'

for (i in 1:length(MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_5$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example5[,1],
          y=example5[,2], 
          z=example5[,3],
          point.col=colors,
          pch = 21,
          sphere.size=1.5,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
view3d(userMatrix = um)

rgl.postscript("./TransitionGroup.3d.3.pdf", "pdf", drawText = FALSE)

label <- findTransitionPath(MarkovObject = MarkovHC_5,
                            level = level,
                            basinA = 1,
                            basinB = 2)
basin[which(label==1)] <- 'transition path\nfrom 1 to 2'
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example5[,1],
          y=example5[,2], 
          z=example5[,3],
          point.col=colors,
          pch = 19,
          sphere.size=1.5,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
view3d(userMatrix = um)
rgl.postscript("./TransitionGroup.3d.4.pdf", "pdf", drawText = FALSE)
```


```R
save.image('./example5.RData')
```


```R
#save plot
pdf('./3D1.basins.pdf')
scatter3D(x=example5[,1],
          y=example5[,2], 
          z=example5[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#08519c",0.5),  alpha("#006d2c",0.5)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D1.path1.pdf')
scatter3D(x=example5[,1],
          y=example5[,2], 
          z=example5[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#08519c",0.5),  alpha("#006d2c",0.5), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D1.path2.pdf')
scatter3D(x=example5[,1],
          y=example5[,2], 
          z=example5[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#08519c",0.5),  alpha("#006d2c",0.5), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D1.path1and2.pdf')
scatter3D(x=example5[,1],
          y=example5[,2], 
          z=example5[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#08519c",0.5),  alpha("#006d2c",0.5), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D1.criticalPoints.pdf')
scatter3D(x=example5[,1],
          y=example5[,2], 
          z=example5[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#08519c",0.5),  alpha("#006d2c",0.5), alpha("#54278f",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()
```
