# simulation data(Fig.2 F-H)


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
library(laGP)
library(ContourFunctions)
```


```R
#theme for Figures
mytheme <-  theme(panel.grid.major =element_blank(),
                  panel.grid.minor = element_blank(),
                  panel.background = element_blank(),
                  axis.line = element_line(size = 1,
                                           colour = "black"),
                  axis.title.x=element_text(size=20,
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

# generate simulation data


```R
#star
example6 <- mvrnorm(n=200, mu=c(10,10), Sigma=matrix(c(5,4.5,4.5,5),2,2))%>%as.data.frame()
example6 <- rbind(example6,
                   mvrnorm(n=200, mu=c(10,10), Sigma=matrix(c(5,-4.5,-4.5,5),2,2)))
example6 <- rbind(example6,
                  mvrnorm(n=200, mu=c(10,10), Sigma=matrix(c(1,0,0,1),2,2)))
#top right Gaussian 
example6 <- rbind(example6,
                   mvrnorm(n=200, mu=c(25,25), Sigma=matrix(c(2,0,0,2),2,2)))
example6 <- rbind(example6,
                  mvrnorm(n=200, mu=c(25,25), Sigma=matrix(c(0.5,0,0,0.5),2,2)))
#low right Gaussian
example6 <- rbind(example6,
                   mvrnorm(n=30, mu=c(25,10), Sigma=matrix(c(3,0,0,3),2,2)))
example6 <- rbind(example6,
                  mvrnorm(n=50, mu=c(25,10), Sigma=matrix(c(0.1,0,0,0.1),2,2)))
#top left Gaussian
example6 <- rbind(example6,
                   mvrnorm(n=300, mu=c(7.5,25), Sigma=matrix(c(1,0,0,1),2,2)))
example6 <- rbind(example6,
                  mvrnorm(n=300, mu=c(13,25), Sigma=matrix(c(1,0,0,1),2,2)))
example6 <- rbind(example6,
                  mvrnorm(n=300, mu=c(10,21), Sigma=matrix(c(1,0,0,1),2,2)))
example6 <- rbind(example6,
                   cbind(runif(n= 20 ,min = 14, max = 18),runif(n= 20 ,min = 25, max = 25)))
example6 <- rbind(example6,
                  cbind(runif(n= 20 ,min = 19, max = 23),runif(n= 20 ,min = 25, max = 25)))
plot(example6)
example6 <- example6[example6[,1]>0,]
```


```R
pdf(file = './ContourAndDensity.pdf', width = 3.5, height = 3.5)
p1 + stat_density_2d(aes(fill = stat(level)), geom = "polygon",
                    h=3,
                    bins = 50,contour=TRUE,alpha=0.001
                    )+scale_fill_gradientn(colours=c(alpha("#252525",0.001)))+mytheme+ theme(legend.position = "none")
dev.off()
```

# Run MarkovHC


```R
MarkovHC_6 <- MarkovHC(origin_matrix=t(example6),
                        transformtype="none",
                        KNN=50,
                        basecluster="kmeans",
                        dobasecluster=TRUE,
                        baseclusternum=400,
                        emphasizedistance=1,
                        weightDist=2,
                        weightDens=0.5,
                        cutpoint=0.01,
                        showprocess=FALSE,
                        bn=2,
                        minBasinSize=0.2,
                        noiseBasinSize=10)
```

# fetch the label of each point


```R
labels <-  fetchLabels(MarkovObject=MarkovHC_6,
                       MarkovLevels=1:length(MarkovHC_6$hierarchicalStructure))

basin <- labels[,23]

for(i in 1:length(basin)){
    basin[i] <- str_split( basin[i], '\\+')[[1]][1]
}

basin <- paste('basin',basin,sep = '')

level <- 23
for (i in 1:length(MarkovHC_6$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[MarkovHC_6$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}

label <- findTransitionPoints(MarkovObject = MarkovHC_6,
                               level = level,
                               basinA = 4,
                               basinB = 1)
basin[which(label==1)] <- 'transition points/nbetween 1 and 4'
dataframe.example6 <- as.data.frame(example6)
dataframe.example6$labels <- factor(basin)
layout <- as.data.frame(dataframe.example6)
```


```R
attractors_layout <- subset(layout, layout$labels=='attractors')
layout <- layout[-which(layout$labels=='attractors'),]
```


```R
ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors_layout, size=1, shape=21, aes(x=V1, y=V2, fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("level 23")+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c("attractors"=alpha("#ce1256",0.8),
               "basin1"=alpha('#8c510a',0.5),
               "basin2"=alpha('#006d2c',0.5),
               "basin3"=alpha('#08519c',0.5),
               "basin4"=alpha('#000000',0.5),
               "basin5"=alpha('#662506',0.5),
               "basin6"=alpha('#993404',0.5)
    ),
    breaks = c("attractors",
               "basin1",
               "basin2",
               "basin3",
               "basin4",
               "basin5",
               "basin6"))+ guides(fill=FALSE)
```


```R
attractors_layout <- subset(layout, layout$labels=='attractors')
layout <- layout[-which(layout$labels=='attractors'),]
# with legend
ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors_layout, size=1, shape=21, aes(x=V1, y=V2, fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("level 23")+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c("attractors"=alpha("#ce1256",0.8),
               "basin1"=alpha('#8c510a',0.5),
               "basin2"=alpha('#006d2c',0.5),
               "basin3"=alpha('#08519c',0.5),
               "basin4"=alpha('#000000',0.5),
               "basin5"=alpha('#662506',0.5),
               "basin6"=alpha('#993404',0.5)
    ),
    breaks = c("attractors",
               "basin1",
               "basin2",
               "basin3",
               "basin4",
               "basin5",
               "basin6"))
```


```R
attractors_layout <- subset(layout, layout$labels=='attractors')
layout <- layout[-which(layout$labels=='attractors'),]

ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors_layout, size=1, shape=21, aes(x=V1, y=V2, fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("level 24")+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c("attractors"=alpha("#ce1256",0.8),
               #"basin"=alpha('#8c510a',0.5),
               "basin1"=alpha('#006d2c',0.5),
               "basin2"=alpha('#08519c',0.5),
               "basin3"=alpha('#000000',0.5),
               "basin4"=alpha('#662506',0.5),
               "basin5"=alpha('#993404',0.5)
    ),
    breaks = c("attractors",
               "basin1",
               "basin2",
               "basin3",
               "basin4",
               "basin5"))+ guides(fill=FALSE)
```


```R
attractors_layout <- subset(layout, layout$labels=='attractors')
layout <- layout[-which(layout$labels=='attractors'),]

ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors_layout, size=1, shape=21, aes(x=V1, y=V2, fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("level 25")+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c("attractors"=alpha("#ce1256",0.8),
               "transition points/nbetween 1 and 4"=alpha('#54278f',1),
               "basin1"=alpha('#006d2c',0.5),
               "basin2"=alpha('#08519c',0.5),
               "basin3"=alpha('#000000',0.5),
               "basin4"=alpha('#662506',0.5)
    ),
    breaks = c("attractors",
               "basin1",
               "basin2",
               "basin3",
               "basin4"))+ guides(fill=FALSE)
```


```R
attractors_layout <- subset(layout, layout$labels=='attractors')
layout <- layout[-which(layout$labels=='attractors'),]

ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors_layout, size=1, shape=21, aes(x=V1, y=V2, fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("level 26")+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c("attractors"=alpha("#ce1256",0.8),
               "transition points/nbetween 1 and 4"=alpha('#54278f',1),
               "basin1"=alpha('#006d2c',0.5),
               "basin2"=alpha('#08519c',0.5),
               "basin3"=alpha('#000000',0.5))
    ,
    breaks = c("attractors",
               "basin1",
               "basin2",
               "basin3"))+ guides(fill=FALSE)
```


```R
centrality_scores <- MarkovHC_6$midResults$centrality_scores

ggplot(data=layout, mapping =  aes(x=V1, y=V2, color=centrality_scores)) +
  geom_point(size=0.5, shape=19)+
  scale_color_gradient2(midpoint=30, 
                        low="#d9d9d9", mid="#969696",high="#000000", space ="Lab" )+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+
  xlab("x") + ylab("y") 
```

# plot the tree


```R
plotHierarchicalStructure(MarkovObject=MarkovHC_6,
                          MarkovLevels=23:26,
                          colorVector=c(alpha('#8c510a',0.5),alpha('#006d2c',0.5),alpha('#08519c',0.5),alpha('#000000',0.5),
                                        alpha('#662506',0.5),alpha('#993404',0.5)))
```

# run MCL


```R
library(MCL)
plot()
adjacency_matrix <- MarkovHC_6[["midResults"]][["symmetric_KNN_graph"]]
adjacency_matrix[adjacency_matrix>0] <- 1
graph_adj <- graph.adjacency(adjacency_matrix,mode="undirected")

pdf(file = './network.pdf', width = 3.5, height = 3.5)
plot(graph_adj,layout=as.matrix(example6),
     vertex.color="#252525",
     vertex.size=1.5,edge.size=0.1,
     edge.color="#969696",vertex.label=NA)#绘出图像
dev.off()

layout <- example6
mcl_results_e2_i0.1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.1, allow1 = TRUE, 
                           max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.2 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.2, allow1 = TRUE, 
                           max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.3 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.3, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.4 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.4, allow1 = TRUE, 
                           max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.5 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.5, allow1 = TRUE, 
                           max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.6 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.6, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.7 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.7, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.8 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.8, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i0.9 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 0.9, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 1, allow1 = TRUE, 
                           max.iter = 100, ESM = FALSE)
mcl_results_e2_i2 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 2, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e2_i3 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 2, inflation = 3, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e1_i1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 1, inflation = 1, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e3_i1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 3, inflation = 1, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e4_i1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 4, inflation = 1, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
mcl_results_e5_i1 <- mcl(x=adjacency_matrix, addLoops = TRUE, expansion = 5, inflation = 1, allow1 = TRUE, 
                         max.iter = 100, ESM = FALSE)
pdf(file = './mcl_results_e1_i1.pdf', width = 3.5, height = 3.5)
layout$labels <- factor(mcl_results_e1_i1$Cluster)
ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("mcl_results_e1_i1.pdf")+
  xlab("x") + ylab("y")+ guides(fill=FALSE)
dev.off()

pdf(file = './mcl_results_e2_i0.6.pdf', width = 3.5, height = 3.5)
layout$labels <- factor(mcl_results_e2_i0.6$Cluster)
ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("mcl_results_e2_i0.6")+
  xlab("x") + ylab("y")+ guides(fill=FALSE)+
  scale_fill_manual(
    values = c("401"=alpha('#006d2c',0.5),
               "1"=alpha('#08519c',0.5),
               "1001"=alpha('#000000',0.5),
               "1400"=alpha('#ae017e',0.5)
               ),
    breaks = c("401",
               "1",
               "1001",'1400'))
dev.off()

pdf(file = './mcl_results_e2_i2.pdf', width = 3.5, height = 3.5)
layout$labels <- factor(mcl_results_e2_i2$Cluster)
ggplot(data=layout, mapping = aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+ggtitle("mcl_results_e2_i2")+
  xlab("x") + ylab("y")+ guides(fill=FALSE)
dev.off()
```
