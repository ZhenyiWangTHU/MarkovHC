setwd("F:/Projects/MarkovHC/Figures/supplementary/spade")
library(spade)
fcsfile <- SPADE.read.FCS(file = "F:/Projects/MarkovHC/data/Cytof/Supplementary Data 2.fcs")
markers <- fcsfile@parameters@data$name
markers <- markers[3:40]
dir.create(paste(getwd(),"/SPADE",sep=""))
SPADE.driver("F:/Projects/MarkovHC/data/Cytof/Supplementary Data 2.fcs",file_pattern="*.fcs",
             out_dir=paste(getwd(),"/SPADE",sep=""),
             cluster_cols=markers,
             panels=NULL,
             comp=TRUE,
             arcsinh_cofactor=NULL,
             transforms=flowCore::arcsinhTransform(a=0, b=0.2),
             downsampling_target_number=NULL,
             downsampling_exclude_pctile=0.02,
             downsampling_target_pctile=0.05,
             k=200,
             clustering_samples=50000,
             layout=igraph:::layout.kamada.kawai,
             pctile_color=c(0.02,0.98))
layout<-read.table(file.path(paste(getwd(),"/SPADE",sep=""),"layout.table"))
mst<-read.graph(file.path(paste(getwd(),"/SPADE",sep=""),"mst.gml"),format="gml")
SPADE.plot.trees(mst,paste(getwd(),"/SPADE",sep=""),file_pattern="*fcs*Rsave",
                 layout=as.matrix(layout),
                 out_dir=file.path(paste(getwd(),"/SPADE",sep=""),"pdf-SPADE"),size_scale_factor=1.2)

