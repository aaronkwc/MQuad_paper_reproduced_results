library(Seurat)

## first read in raw gene expression matrix downloaded from GEO
tnbc.counts <- read.table('/storage/yhhuang/research/mito/copyKAT/tnbc1/GSM4476486_filtered_UMIcount_TNBC1.txt')

## keep tumor cells only
# use copyKAT results and remove clone 2 (normal cells)
copykat <- read.table('~/data/reproduced_results/droplet_based_data/TNBC1_data/TNBC1_copykat.txt')
keep <- copykat$V1 != 2
tnbc.counts <- tnbc.counts[,keep]

##run standard seurat workflow on the downloaded data
tnbc <- CreateSeuratObject(counts = tnbc.counts)
tnbc <- NormalizeData(object = tnbc)
tnbc <- FindVariableFeatures(object = tnbc)
tnbc <- ScaleData(object = tnbc)
tnbc <- RunPCA(object = tnbc)
ElbowPlot(tnbc)
tnbc <- FindNeighbors(object = tnbc)
tnbc <- FindClusters(object = tnbc)
tnbc <- RunUMAP(object = tnbc, dims = 1:20)
DimPlot(object = tnbc, reduction = "umap")


## read in mquad clones for DE analysis
mquad <- read.csv('~/data/reproduced_results/droplet_based_data/TNBC1_data/clones_df.csv')
clone.id <- mquad$clone_id[mquad$copykat != 2]
tnbc$mquad_id <- clone.id

tnbc <- SetIdent(tnbc, value = tnbc$mquad_id)
DimPlot(object = tnbc, reduction = "umap")

t1.markers <- FindMarkers(tnbc, ident.1 = 1, min.pct = 0.25, test.use = "DESeq2")
head(t1.markers)

## save umap cooordinates and markers to python for visualization with mquad clones
write.csv(tnbc@reductions[["umap"]]@cell.embeddings, file='~/data/reproduced_results/droplet_based_data/TNBC1_data/TNBC1_UMAP.csv')
write.csv(t1.markers, file='~/data/reproduced_results/droplet_based_data/TNBC1_data/TNBC1_t1_markers.csv')
