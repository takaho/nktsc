# Seurat macro
# Rscript seurat.R [filtered_sparse_matrix_directory] [output directory]
library(dplyr)
library(Seurat)

library(reticulate)
use_python('/Users/takaho/.pyenv/shims/python')
srcdir<-"."
outdir<-"seurat_out"
args<-commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
    srcdir<-args[1]
}
if (length(args) > 1) {
    outdir<-args[2]
}
if (!dir.exists(outdir)) {
    dir.create(outdir, recursive=TRUE)
}

output <- file.path(outdir, "seurat.data.rds")
umap_file <- file.path(outdir, "seurat.umap.tsv")
normalizedfile <- file.path(outdir, "seurat.normalized.tsv")
countfile <- file.path(outdir, "seurat.count.tsv")
filename_pca <- file.path(outdir, "seurat.pca.tsv")
filename_cluster <- file.path(outdir, "seurat.clusters.tsv")
filename_mito = file.path(outdir, "seurat.mt.tsv")

if (length(args) > 2) {
    output<-args[2]
}
if (length(args) > 3) {
    umap_file<-args[4]
}

print("READ DATA")
sc.data <- Read10X(data.dir = srcdir) 

# Filtering using mitchondria genes were not applied in the study to identify clusters of dead cells 
#sc.data <- Read10X(data.dir = srcdir, gene.column=1) # if only 1 column in features.tsv
#dim(sc.data)

# Initialize the Seurat object with the raw (non-normalized data).
sc <- CreateSeuratObject(counts = sc.data, project=outdir, min.cells = 5, min.features = 200)
sc

# print("Find mt-genes")
# sc[["percent.mt"]] <- PercentageFeatureSet(sc, pattern = "^mt-") # if human use MT
#write.table(sc[["percent.mt"]], file=filename_mito, sep='\t')

print("Normalize data wigh LogNormalize method")
sc <- subset(sc, subset = nFeature_RNA > 200 & nFeature_RNA < 20000)
write.table(sc[["RNA"]]@data, countfile, sep="\t")
sc<- NormalizeData(sc, normalization.method = "LogNormalize", scale.factor = 10000)
sc

write.table(sc[["RNA"]]@data, normalizedfile, sep='\t')

sc <- FindVariableFeatures(sc, selection.method = "vst", nfeatures = 2000)

# Identify the 10 most highly variable genes
print("TOP 10")
top10 <- head(VariableFeatures(sc), 10)
print(length(VariableFeatures(sc)))
top10

print("PCA")
all.genes <- rownames(sc)
sc <- ScaleData(sc, features = all.genes)

sc <- RunPCA(sc, features = VariableFeatures(object = sc), npcs=50)

VizDimLoadings(sc, dims = 1:2, reduction = "pca")

write.table(attributes(sc[['pca']])$cell.embeddings, file=filename_pca, sep='\t')

print("clusters")
sc <- FindNeighbors(sc, dims = 1:20)
sc <- FindClusters(sc, resolution = 0.6)

write.table(Idents(sc), filename_cluster, sep="\t")

print("UMAP")
sc <- RunUMAP(sc, n.components=3, dims = 1:10 , umap.method='umap-learn', metric='correlation')

write.table(attributes(sc$umap)$cell.embedding, file=umap_file, sep="\t")
