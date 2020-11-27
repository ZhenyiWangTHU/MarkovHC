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

# generate simulation data


```R
#two basins have high density
example4 <- mvrnorm(n=1000, mu=c(10,10,10), Sigma=matrix(c(10,0,0,0,10,0,0,0,10),3,3))%>%as.data.frame()

example4 <- rbind(example4,
                  mvrnorm(n=1000, mu=c(50,10,10), Sigma=matrix(c(10,0,0,0,10,0,0,0,10),3,3)))

example4 <- rbind(example4,
                  cbind(runif(n= 50 ,min = 10, max = 27),
                        runif(n= 50 ,min = 10, max = 10),
                        runif(n= 50 ,min = 10, max = 10)))

example4 <- rbind(example4,
                  cbind(runif(n= 50 ,min = 32, max = 50),
                        runif(n= 50 ,min = 10, max = 10),
                        runif(n= 50 ,min = 10, max = 10)))
#usethis::use_data(example4,overwrite = FALSE)

scatter3d(x=example4[,1],
          y=example4[,2], 
          z=example4[,3],
          point.col='black',
          pch = 21,
          sphere.size=2,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE)
```

# run MarkovHC


```R
MarkovHC_4 <- MarkovHC(origin_matrix=t(example4),
                       transformtype="none",
                       KNN=30,
                       basecluster="kmeans",
                       dobasecluster=TRUE,
                       baseclusternum=100,
                       emphasizedistance=1,
                       weightDist=2,
                       weightDens=0.5,
                       cutpoint=0.01,
                       showprocess=FALSE,
                       bn=2,
                       minBasinSize=0.2,
                       noiseBasinSize=10)
```


```R
labels <-  fetchLabels(MarkovObject=MarkovHC_4,
                       MarkovLevels=1:length(MarkovHC_4$hierarchicalStructure))

basin <- labels[,18]

for(i in 1:length(basin)){
    basin[i] <- str_split( basin[i], '\\+')[[1]][1]
}

basin <- paste('basin',basin,sep = '')
level <- 18
for (i in 1:length(MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
```


```R
#find the transition points
label <- findTransitionPoints(MarkovObject = MarkovHC_4,
                              level = level,
                              basinA = 1,
                              basinB = 2)
basin[which(label==1)] <- 'transition points\nbetween 1 and 2'

colors <- c("#e41a1c","#377eb8","#4daf4a","#984ea3")
colors <- colors[factor(basin)]
scatter3d(x=example4[,1],
          y=example4[,2], 
          z=example4[,3],
          point.col=colors,
          pch = 21,
          sphere.size=2,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
rgl.postscript("./TransitionGroup.3d.1.pdf", "pdf", drawText = FALSE)

#find the transition paths
level <- 18
label <- findTransitionPath(MarkovObject = MarkovHC_4,
                            level = level,
                            basinA = 1,
                            basinB = 2)
basin[which(label==1)] <- 'transition path\nfrom 1 to 2'

for (i in 1:length(MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example4[,1],
          y=example4[,2], 
          z=example4[,3],
          point.col=colors,
          pch = 21,
          sphere.size=2,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
rgl.postscript("./TransitionGroup.3d.2.pdf", "pdf", drawText = FALSE)

level <- 18
label <- findTransitionPath(MarkovObject = MarkovHC_4,
                            level = level,
                            basinA = 2,
                            basinB = 1)
basin[which(label==1)] <- 'transition path\nfrom 2 to 1'

for (i in 1:length(MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example4[,1],
          y=example4[,2], 
          z=example4[,3],
          point.col=colors,
          pch = 21,
          sphere.size=2,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
rgl.postscript("./TransitionGroup.3d.3.pdf", "pdf", drawText = FALSE)

#transition path
level <- 18
label <- findTransitionPath(MarkovObject = MarkovHC_4,
                            level = level,
                            basinA = 1,
                            basinB = 2)
basin[which(label==1)] <- 'transition path\nbetween 1 and 2'
label <- findTransitionPath(MarkovObject = MarkovHC_4,
                            level = level,
                            basinA = 2,
                            basinB = 1)
basin[which(label==1)] <- 'transition path\nbetween 1 and 2'

for (i in 1:length(MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_4$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
colors <- c("#e41a1c","#377eb8","#4daf4a","#ff7f00")
colors <- colors[factor(basin)]
scatter3d(x=example4[,1],
          y=example4[,2], 
          z=example4[,3],
          point.col=colors,
          pch = 21,
          sphere.size=2,
          surface=FALSE,
          xlab = "X", ylab = "Y",
          zlab = "Z", axis.scales = FALSE,
          color = colors)
rgl.postscript("./TransitionGroup.pathbetween1and2.pdf", "pdf", drawText = FALSE)
```


```R
save.image('./example4.RData')
```


```R
#使用没有景深的三维图
pdf('./3D2.basins.pdf')
scatter3D(x=example4[,1],
          y=example4[,2], 
          z=example4[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#006d2c",0.2),  alpha("#08519c",0.2)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D点2.path1.pdf')
scatter3D(x=example4[,1],
          y=example4[,2], 
          z=example4[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#006d2c",0.2),  alpha("#08519c",0.2), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D2.path2.pdf')
scatter3D(x=example4[,1],
          y=example4[,2], 
          z=example4[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#006d2c",0.2),  alpha("#08519c",0.2), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D2.path1and2.pdf')
scatter3D(x=example4[,1],
          y=example4[,2], 
          z=example4[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#006d2c",0.2),  alpha("#08519c",0.2), alpha("#feb24c",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()

pdf('./3D2.criticalPoints.pdf')
scatter3D(x=example4[,1],
          y=example4[,2], 
          z=example4[,3], 
          colvar = as.numeric(as.factor(basin)), 
          col = c(alpha("#ce1256",0.7),alpha("#006d2c",0.2),  alpha("#08519c",0.2), alpha("#54278f",1)),
          pch = 19, cex = 1.5,theta = 40, phi = 40)
dev.off()
```
