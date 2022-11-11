## QC and filtering 

## Prepare data for quality metrics

load("merge_obj.RData")


View(merge_obj@meta.data)

# create a sample column
merge_obj$sample <- rownames(merge_obj@meta.data)

# split sample column
merge_obj@meta.data <- separate(merge_obj@meta.data, col = 'sample', into = c('Sample_no','Barcode'), 
                                sep = '_')
## Read phenotype data to R
phen=read.delim2('phendata.txt', sep = ',', header = T)

phen=phen %>%
  filter(Sample.Name %in% c('GSM5676879','GSM5676880','GSM5676883', 'GSM5676884'))

## integrate phendata with merge_obj
for (y in 1:nrow(merge_obj@meta.data)){
  z=merge_obj@meta.data$Sample_no[y]
  merge_obj@meta.data$sex[y] = phen$sex[phen$Sample.Name==z]
  merge_obj@meta.data$source_name[y] = phen$source_name[phen$Sample.Name==z]
  merge_obj@meta.data$phenotype[y] = phen$Phenotype[phen$Sample.Name==z]
}



### Perform quality assesment of the data
# calculate mitochondrial percentage
merge_obj$mitoPercent <- PercentageFeatureSet(merge_obj, pattern='^MT-')
merge_obj$riboPercent <- PercentageFeatureSet(merge_obj, pattern="^RP[LS]")

## Get a df from the slot meta data from the seurat object 

merge_obj@meta.data -> meta.data


## Add some colunms to the df 

# Add number of genes per UMI for each cell to metadata
meta.data$log10GenesPerUMI <- log10(meta.data$nFeature_RNA) / log10(meta.data$nCount_RNA)

# Compute percent mito ratio
#meta.data$mitoRatio <- PercentageFeatureSet(object = merge_obj, pattern = "^MT-")
meta.data$mitoRatio <- meta.data$mitoPercent / 100

merge_obj@meta.data <- meta.data 



# Clone obj to perform some manipulation 
merge_obj -> data.nomalat

#
#*****23762 genes and 30587 cells in 4 samples in its sparse matrix Only nonezeros values are stored
#Compute the percentage of the largest gene 
apply(data.nomalat@assays$RNA@counts,2,max) -> data.nomalat$largest_count

apply(data.nomalat@assays$RNA@counts,2,which.max)-> data.nomalat$largest_index

rownames(data.nomalat)[data.nomalat$largest_index] -> data.nomalat$largest_gene

100 * data.nomalat$largest_count / data.nomalat$nCount_RNA -> data.nomalat$percent.Largest.Gene

data.nomalat$largest_gene -> merge_obj$largest_gene
data.nomalat$percent.Largest.Gene -> merge_obj$percent.Largest.Gene
rm(data.nomalat)

merge_obj@meta.data -> meta.data 

# Rename columns
meta.data <- meta.data %>%
  dplyr::rename(seq_folder = Sample_no,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)


# Add metadata back to Seurat object
merge_obj@meta.data <- meta.data


# Create .RData object to load at any time
save(merge_obj, 
     file=paste0(paths,'/','merged_obj.RData'))

#################################################################
#load the merged_obj
library(ggplot2)

load('merged_obj.RData')

merge_obj@meta.data -> meta.data
## Exploring the quality of the data ..........................................


# Visualize the number of cell counts per sample
meta.data %>% 
  ggplot(aes(x=phenotype, fill=phenotype)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

# Visualize the number UMIs/transcripts per cell
meta.data %>% 
  ggplot(aes(color=phenotype, x=nUMI, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 23000)



# Visualize the distribution of genes detected per cell via histogram
meta.data %>% 
  ggplot(aes(color=phenotype, x=nGene, fill= phenotype)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300)


# filtering# Visualize the distribution of genes detected per cell via boxplot
meta.data %>% 
  ggplot(aes(x=phenotype, y=log10(nGene), fill=phenotype)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

meta.data %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoPercent)) + 
  geom_point() + 
  scale_colour_gradient(low = "grey", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept =1000) +
  geom_hline(yintercept = 500) +
  facet_wrap(~phenotype)


meta.data %>% 
  ggplot(aes(color=phenotype, x=mitoPercent, fill=phenotype)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic()


# Visualize the overall complexity of the gene expression by visualizing the genes detected per UMI
meta.data %>%
  ggplot(aes(x=log10GenesPerUMI, color = phenotype, fill=phenotype)) +
  geom_density(alpha = 0.2) +
  theme_classic()

VlnPlot(merge_obj, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))


## Cell level filtering 
merge_filtered <- subset(merge_obj, subset = nUMI > 1000 &
                           nGene > 500 &
                           mitoPercent > 0.5 & 
                           log10GenesPerUMI > 0.85)

# View quality after 
VlnPlot(merge_filtered, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))


merge_filtered

#merged_seurat



## gene level filtering 


# Output a logical vector for every gene on whether the more than zero counts per cell
# Extract counts
counts <- GetAssayData(object = merge_filtered, slot = "counts")

# Output a logical vector for every gene on whether the more than zero counts per cell
nonzero <- counts > 0

# Sums all TRUE values and returns TRUE if more than 10 TRUE values per gene
keep_genes <- Matrix::rowSums(nonzero) >= 10

# Only keeping those genes expressed in more than 10 cells
filtered_counts <- counts[keep_genes, ]

# Reassign to filtered Seurat object
filtered_AMDP <- CreateSeuratObject(filtered_counts, meta.data = merge_filtered@meta.data)



# Visualiz 
VlnPlot(filtered_AMDP, features=c("nUMI","mitoPercent", "riboPercent","percent.Largest.Gene"))

#filtered_AMDP@meta.data->meta

#meta$source_name[grep('Chorod', meta$source_name)]='Choroid'

#merge_obj@meta.data <- meta.data

#filtered_AMDP@meta.data<-meta


# perform standard workflow steps to figure out if we see any batch effects --------

## Normalize the data 
merged_seurat_filtered <- NormalizeData(object = filtered_AMDP)
##Feature selection
merged_seurat_filtered <- FindVariableFeatures(object = merged_seurat_filtered)
## Scaling data
merged_seurat_filtered <- ScaleData(object = merged_seurat_filtered)
## Dimensional Reduction
merged_seurat_filtered <- RunPCA(object = merged_seurat_filtered)

ElbowPlot(merged_seurat_filtered)
##Cell Sub population Identification using unsupervised clustering methods [!prior information based]
## Techniques for unsupervised clustering 
## (i) k-means; (ii) hierarchical clustering; (iii) density-based clustering; and (iv) graph-based clustering 
## 
## Other cluster methods are 
## single-cell consensus clustering (SC3) (Kiselev et al., 2017) 
## shared nearest neighbor (SNN) clustering algorithm (Satija et al., 2015) ....Seurat[findNeighbors]
##SC3 is an unsupervised approach that combines multiple clustering approaches, 
## which has a high accuracy and robustness in single-cell clustering.
merged_seurat_filtered <- FindNeighbors(object = merged_seurat_filtered, dims = 1:20)
merged_seurat_filtered <- FindClusters(object = merged_seurat_filtered)
merged_seurat_filtered <- RunUMAP(object = merged_seurat_filtered, dims = 1:20)


# plot
p1 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'phenotype')
p2 <- DimPlot(merged_seurat_filtered, reduction = 'umap', group.by = 'source_name',
              cols = c('red','green'))

grid.arrange(p1, p2, ncol = 2, nrow = 2)
#pca plot 
t1 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'phenotype')
t2 <- DimPlot(merged_seurat_filtered, reduction = 'pca', group.by = 'source_name',
              cols = c('red','green'))

## Get the list of most highly expressed genes
apply(merged_seurat_filtered@assays$RNA@data,1,mean) -> gene.expression

sort(gene.expression, decreasing = TRUE) -> gene.expression

head(gene.expression, n=10)

save(filtered_AMDP, file = 'AMD_Project.RData')
