##################################################################
## Check the impact of cell cycle on the expression score 
##Predict the cell cycle of each cell 
## Get a new obj that allows for examine the cell cycle phase 

CellCycleScoring(seurat.integrated, 
                 s.features = cc.genes.updated.2019$s.genes, 
                 g2m.features = cc.genes.updated.2019$g2m.genes, 
                 set.ident = TRUE) -> cell.data

View(cell.data@meta.data)

as_tibble(cell.data[[]]) %>%
  ggplot(aes(Phase)) + geom_bar()

## Another form of visual 
as_tibble(cell.data[[]]) %>%
  ggplot(aes(x=S.Score, y=G2M.Score, color=Phase)) + 
  geom_point() +
  coord_cartesian(xlim=c(-0.15,0.15), ylim=c(-0.15,0.15))


### Most of the cells in are in G1 phase which means they in cell growth phase. 

## Getting  the top 500 most varied genes 
FindVariableFeatures(
  seurat.integrated, 
  selection.method = "vst", 
  nfeatures=500
) -> var.data


as_tibble(HVFInfo(var.data),rownames = "Gene") -> variance.data

variance.data %>% 
  mutate(hypervariable=Gene %in% VariableFeatures(var.data)
  ) -> variance.data

head(variance.data, n=10)

## Compare in graph
variance.data %>% 
  ggplot(aes(log(mean),log(variance),color=hypervariable)) + 
  geom_point() + 
  scale_color_manual(values=c("black","red"))


####**************************************DIMENSION REDUCTION ****************#######
## PCA to reduce the dimension 
RunPCA(seurat.integrated,features=VariableFeatures(seurat.integrated)) -> PCA.data

DimPlot(PCA.data, reduction = 'pca')

DimPlot(PCA.data,reduction="pca", group.by = "phenotype", 
        label = F, label.size = 3)

ElbowPlot(PCA.data)

DimHeatmap(PCA.data,dims=1:6, cells=500)


#### Performimg tSNE
8482 -> saved.seed
set.seed(saved.seed)

save(PCA.data, file = "PCA.RData")


RunTSNE(
  PCA.data,
  dims=1:15,
  seed.use = saved.seed, 
  perplexity=100
) -> tdata

DimPlot(tdata,reduction = "tsne", pt.size = 1)+
  ggtitle("tSNE with Perplexity 100")


#save

save(tdata, file = "tdata.RData")
