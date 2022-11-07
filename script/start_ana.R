#Library----
library(Seurat)
library(dplyr)
library(khroma)

#Import data----
Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm.RDS")
# Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_log.RDS")

my_palette = c("VV-B8R"= "#DDCC77","Ly49p-CD8ab"= "#117733","VV-Ly49p-CD8aa"= "#882255",
  "Ly49p-CD8aa"= "#88CCEE","VV-Ly49p-CD8ab"= "#999933","Naive"= "#AA4499",
  "Ly49n"= "#332288")

#PCA----
Singlet_norm <- RunPCA(Singlet_norm, features = VariableFeatures(object = Singlet_norm), assay = "SCT")
print(Singlet_norm[["pca"]], dims = 1:5, nfeatures = 5)
VizDimLoadings(Singlet_norm, dims = 1:4, reduction = "pca", nfeatures = 100)
VizDimLoadings(Singlet_norm, dims = 1, reduction = "pca", nfeatures = 40)
ElbowPlot(Singlet_norm, ndims = 50)
DimPlot(Singlet_norm, reduction = "pca", group.by = "HTO_maxID")
DimHeatmap(Singlet_norm, dims = 1:6, cells = 500, balanced = TRUE, nfeatures = 20)
DimHeatmap(Singlet_norm, dims = 7:12, cells = 500, balanced = TRUE, nfeatures = 20)

#UMAP----
Singlet_norm <- RunUMAP(Singlet_norm, dims = 1:13, assay = "SCT")
DimPlot(Singlet_norm, reduction = "umap",group.by = "HTO_maxID", pt.size = 1, cols = my_palette)
