##****************Cell type identification with scType**********************************

## Cell type assignment 
## Load the scType functions 
# load gene set preparation function

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/gene_sets_prepare.R")
# load cell type annotation function
source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_score_.R")

# To prepare the gene set lets import the db

# DB file
db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
tissue = "Eye" # e.g. Immune system,Pancreas,Liver,Eye,Kidney,Brain,Lung,Adrenal,Heart,Intestine,Muscle,Placenta,Spleen,Stomach,Thymus 

# prepare gene sets
gs_list = gene_sets_prepare(db_, tissue)


# Lets assign the cell  types to each cluster 

es.max = sctype_score(scRNAseqData = tfdata@assays$integrated@scale.data, scaled = TRUE,
                      gs = gs_list$gs_positive, gs2 = gs_list$gs_negative)

# Merge by cluster

cL_results = do.call("rbind", lapply(unique(tfdata@meta.data$seurat_clusters),
                                     function(cl){
                                       es.max.cl = sort(rowSums(es.max[,row.names(tfdata@meta.data[tfdata@meta.data$seurat_clusters==cl, ])]),
                                                        decreasing = !0)
                                       head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(tfdata@meta.data$seurat_clusters==cl)),
                                            10)
                                       
                                     }))

sctype_scores = cL_results %>% 
  group_by(cluster) %>% 
  top_n(n=1,wt =scores)

# set low-confident (low ScType score) clusters to "unknown"


sctype_scores$type[as.numeric(as.character(sctype_scores$scores)) < sctype_scores$ncells /
                     4] = "Unknown"

View(sctype_scores[, 1:3])


# Lets visualize the assigned cell 

tfdata@meta.data$Cell_Identity = " "
for(j in unique(sctype_scores$cluster)){
  cl_type = sctype_scores[sctype_scores$cluster==j,];
  tfdata@meta.data$Cell_Identity[tfdata@meta.data$seurat_clusters == j] = as.character(cl_type$type[1])
}

DimPlot(tfdata, reduction = "tsne", label = T, repel = T, group.by = "Cell_Identity")

DimPlot(tfdata, reduction = "umap", label = T, repel = T, group.by = c("Cell_Identity", 
                                                                       "phenotype"))

DimPlot(tfdata, reduction = "umap", label = T, repel = T, group.by = "Cell_Identity")

## Visualise bubble plot showing all the cell types that were considered by scType 

# load libraries
lapply(c("ggraph","igraph","tidyverse", "data.tree"), library, character.only = T)

# prepare edges
cL_results=cL_results[order(cL_results$cluster),]; edges = cL_results; edges$type = paste0(edges$type,"_",edges$cluster); edges$cluster = paste0("cluster ", edges$cluster); edges = edges[,c("cluster", "type")]; colnames(edges) = c("from", "to"); rownames(edges) <- NULL

# prepare nodes
nodes_lvl1 = sctype_scores[,c("cluster", "ncells")]; 
nodes_lvl1$cluster = paste0("cluster ", nodes_lvl1$cluster); 
nodes_lvl1$Colour = "#f1f1ef"; nodes_lvl1$ord = 1; 
nodes_lvl1$realname = nodes_lvl1$cluster; 
nodes_lvl1 = as.data.frame(nodes_lvl1); 
nodes_lvl2 = c();
ccolss= c("#5f75ae","#92bbb8","#64a841","#e5486e","#de8e06","#eccf5a","#b5aa0f","#e4b680","#7ba39d","#b15928","#ffff99", "#6a3d9a","#cab2d6","#ff7f00","#fdbf6f","#e31a1c","#fb9a99","#33a02c","#b2df8a","#1f78b4","#a6cee3")
for (i in 1:length(unique(cL_results$cluster))){
  dt_tmp = cL_results[cL_results$cluster == unique(cL_results$cluster)[i], ]; nodes_lvl2 = rbind(nodes_lvl2, data.frame(cluster = paste0(dt_tmp$type,"_",dt_tmp$cluster), ncells = dt_tmp$scores, Colour = ccolss[i], ord = 2, realname = dt_tmp$type))
}
nodes = rbind(nodes_lvl1, nodes_lvl2); nodes$ncells[nodes$ncells<1] = 1;
files_db = openxlsx::read.xlsx(db_)[,c("cellName","shortName")]; files_db = unique(files_db); nodes = merge(nodes, files_db, all.x = T, all.y = F, by.x = "realname", by.y = "cellName", sort = F)
nodes$shortName[is.na(nodes$shortName)] = nodes$realname[is.na(nodes$shortName)]; nodes = nodes[,c("cluster", "ncells", "Colour", "ord", "shortName", "realname")]

library(igraph)
# remove duplicates from cluster colume
nodess = nodes[!duplicated(nodes$cluster),]
# Create a graph object from the data frame
mygraph <- graph_from_data_frame(edges, vertices=nodess)


plot.igraph(mygraph)


library(ggraph)
# Make the graph
gggr<- ggraph(mygraph, layout = 'circlepack', weight=I(ncells)) + 
  geom_node_circle(aes(filter=ord==1,fill=I("#F5F5F5"), colour=I("#D3D3D3")), alpha=0.9) + geom_node_circle(aes(filter=ord==2,fill=I(Colour), colour=I("#D3D3D3")), alpha=0.9) +
  theme_void() + geom_node_text(aes(filter=ord==2, label=shortName, colour=I("#ffffff"), fill="white", repel = !1, parse = T, size = I(log(ncells,25)*1.5)))+ geom_node_label(aes(filter=ord==1,  label=shortName, colour=I("#000000"), size = I(3), fill="white", parse = T), repel = !0, segment.linetype="dotted")


# Define multiplot function 
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}
library(ggplot2)
multiplot(DimPlot(tfdata, reduction = "umap", label = TRUE, repel = TRUE, cols = ccolss), gggr, cols = 2)
