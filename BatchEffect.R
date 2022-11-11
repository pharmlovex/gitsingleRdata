


load("AMD_Project.RData")

###**************NORMALISATION & BATCH EFFECT CORRECTION*****************

# perform integration to correct for batch effects ------
#obj.list <- SplitObject(merged_seurat_filtered, split.by = 'phenotype')

obj.list <- SplitObject(filtered_AMDP, split.by = 'phenotype')

for(i in 1:length(obj.list)){
  obj.list[[i]] <- NormalizeData(object = obj.list[[i]])
  obj.list[[i]] <- FindVariableFeatures(object = obj.list[[i]])
}
# select integration features
features <- SelectIntegrationFeatures(object.list = obj.list)

# find integration anchors (CCA)
anchors <- FindIntegrationAnchors(object.list = obj.list,
                                  anchor.features = features, reduction = "cca")


save(anchors, file="anchors.RData")

load("anchors.RData")
# integrate data
seurat.integrated <- IntegrateData(anchorset = anchors)


# Scale data, run PCA and UMAP and visualize integrated data
seurat.integrated <- ScaleData(object = seurat.integrated)
seurat.integrated <- RunPCA(object = seurat.integrated)
seurat.integrated <- RunUMAP(object = seurat.integrated, dims = 1:50)


p3 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'phenotype')
p4 <- DimPlot(seurat.integrated, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))
p3

grid.arrange(p1, p2, p3, p4, ncol = 2, nrow = 2)


t3=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'phenotype')
t4=DimPlot(seurat.integrated, reduction = 'pca', group.by = 'source_name',
           cols = c('red','green'))

## View all the pca plots
grid.arrange(t1, t3, t2, t4, ncol = 2, nrow = 2)


#save

save(seurat.integrated, file = 'AMDIntregate.RData')

