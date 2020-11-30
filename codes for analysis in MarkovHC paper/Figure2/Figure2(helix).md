```R
#simulate data
library(MarkovHC)
library(dplyr)
library(car)
library(dplyr)
library(MASS)
library(ggplot2)

data("helix_data", package = "MarkovHC")
```


```R
#Figure theme
mytheme <-  theme(panel.grid.major =element_blank(),
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
                                              face = "bold"),
                  legend.background = element_blank(),
                  legend.key=element_blank())
```

# run MarkovHC


```R
Markov_helix <- MarkovHC(origin_matrix=t(helix_data),
                       transformtype="none",
                       KNN=30,
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
                       noiseBasinSize=10)
```


```R
labels <-  fetchLabels(MarkovObject=Markov_helix,
                       MarkovLevels=1:length(Markov_helix$hierarchicalStructure))

basin <- labels[,30]

for(i in 1:length(basin)){
    basin[i] <- str_split( basin[i], '\\+')[[1]][1]
}

basin <- paste('basin',basin,sep = '')
level <- 30
for (i in 1:length(Markov_helix$hierarchicalStructure[[level]]$attractorPoints)) {
  basin[Markov_helix$hierarchicalStructure[[level]]$attractorPoints[[i]]] <- "attractors"
}
```


```R
dataframe.helix_data <- as.data.frame(helix_data)
dataframe.helix_data$labels <- factor(basin)
layout <- as.data.frame(dataframe.helix_data)
```


```R
attractors <- subset(layout, layout$labels=="attractors")
layout <- layout[-which(layout$labels=='attractors'),]

ggplot(data=layout, mapping =  aes(x=V1, y=V2)) +
  geom_point(size=1, shape=21, aes(fill=labels), color=alpha("#525252",0))+
  geom_point(data=attractors,size=1, shape=21, aes(x=V1, y=V2, fill=labels), 
             color=alpha("#525252",0))+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+
  xlab("x") + ylab("y")+
  scale_fill_manual(
    values = c(  "basin1"=alpha('#006d2c',0.5),
                 "basin2"=alpha('#08519c',0.5),
                 "basin3"=alpha('#000000',0.5),
                 "attractors"=alpha("#ce1256",0.8)
               ),
    breaks = c("basin1",
               "basin2",
               "basin3",
               "attractors"
               ))+ guides(fill=FALSE)
```


```R
layout <- as.data.frame(dataframe.helix_data)
centrality_scores <- Markov_helix$midResults$centrality_scores

ggplot(data=layout, mapping =  aes(x=V1, y=V2, color=centrality_scores)) +
  geom_point(size=0.5, shape=19)+
  scale_color_gradient2(midpoint=20, 
                        low="#d9d9d9", mid="#969696",high="#000000", space ="Lab" )+
  xlim(min(layout$V1)-1,max(layout$V1)+1)+
  ylim(min(layout$V2)-1,max(layout$V2)+1)+
  mytheme+
  xlab("x") + ylab("y") #+ theme(legend.position = 'none')
```


```R
save.image('./helix.RData')
```
