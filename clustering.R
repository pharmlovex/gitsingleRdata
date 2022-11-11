



load("tdata.RData")


###Defining Cluster using graph based methods
FindNeighbors(tdata,dims=1:15) -> kdata  ##--------kmeans = 20 
## Distance to the nearest 20 neigbours
kdata@graphs$integrated_snn[1:10,1:10]
##Segment the graph with FindCluster 
FindClusters(kdata,resolution = 0.03) -> fdata


###Cluster metadata 

DimPlot(fdata,reduction="pca",label = TRUE)+
  ggtitle("PC1 vs PC2 with Clusters")


DimPlot(fdata,reduction="pca", dims=c(4,9), label=TRUE)+
  ggtitle("PC4 vs PC9 with Clusters")


#### Performimg tSNE
8482 -> saved.seed
set.seed(saved.seed)


## check tSNE plot 
RunTSNE(
  fdata,
  dims=1:15,
  seed.use = saved.seed, 
  perplexity=100
) -> tfdata

DimPlot(tfdata,reduction="tsne",
        pt.size = 1, label = TRUE, label.size = 7)






## Compute the QC metric for thr clustered data 
VlnPlot(fdata,features="nCount_RNA")

VlnPlot(fdata,features="nFeature_RNA")

VlnPlot(fdata,features= 'mitoPercent', group.by = 'phenotype')

## Cell cycle per cluster 
fdata@meta.data %>%
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Percentage of cell cycle phases per cluster")


 
## Cell cycle per cluster - Normal
fdata@meta.data %>%
  filter(phenotype == "Normal") %>% 
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Normal Eye") -> cellCylNorm
## Cell cycle per cluster - Disease
fdata@meta.data %>%
  filter(phenotype == "Early AMD") %>% 
  group_by(seurat_clusters,Phase) %>%
  count() %>%
  group_by(seurat_clusters) %>%
  mutate(percent=100*n/sum(n)) %>%
  ungroup() %>%
  ggplot(aes(x=seurat_clusters,y=percent, fill=Phase)) +
  geom_col() +
  ggtitle("Early AMD") -> cellcylAMD



cellCylNorm|cellcylAMD


fdata[[]] %>%
  group_by(seurat_clusters, phenotype) %>%
  count() %>%
  arrange(desc(n)) %>%
  group_by(seurat_clusters) %>%
  slice(1:2) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(n))


VlnPlot(fdata,features="MALAT1")

VlnPlot(fdata,features="percent.Largest.Gene", group.by = 'phenotype')

## Which gene is the largest 
fdata[[]] %>%
  filter(largest_gene != 'MALAT1') %>%
  group_by(seurat_clusters, largest_gene) %>%
  count() %>%
  arrange(desc(n)) %>%
  group_by(seurat_clusters) %>%
  slice(1:2) %>%
  ungroup() %>%
  arrange(seurat_clusters, desc(n)) -> data.gene



tfdata@reductions$tsne@cell.embeddings %>%
  as_tibble() %>%
  add_column(seurat_clusters=tfdata$seurat_clusters, largest_gene=tfdata$largest_gene) %>%
  filter(largest_gene %in% largest_genes_to_plot) %>%
  ggplot(aes(x=tSNE_1, y=tSNE_2, colour=seurat_clusters)) +
  geom_point() +
  facet_wrap(vars(largest_gene))



