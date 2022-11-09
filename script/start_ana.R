#Library----
library(Seurat)
library(dplyr)
library(khroma)
library(future) # parallelization in Seurat (NormalizeData() / ScaleData() / JackStraw() / FindMarkers() / FindIntegrationAnchors() / FindClusters() - if clustering over multiple resolutions)
library(clustree)
library(ggplot2)


#Start session----
rm(list = ls())
plan("multiprocess", workers = 9) # activate parallelization

#Import data----
Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm.RDS")
# Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_log.RDS")

my_palette = c("VV-B8R"= "#DDCC77",
               "Ly49p-CD8ab"= "#117733",
               "VV-Ly49p-CD8aa"= "#882255",
               "Ly49p-CD8aa"= "#88CCEE",
               "VV-Ly49p-CD8ab"= "#999933",
               "Naive"= "#AA4499",
               "Ly49n"= "#332288")

#PCA----
Singlet_norm <- RunPCA(Singlet_norm, features = VariableFeatures(object = Singlet_norm), assay = "SCT")
print(Singlet_norm[["pca"]], dims = 1:4, nfeatures = 10)
VizDimLoadings(Singlet_norm, dims = 1:4, reduction = "pca", nfeatures = 20)
ElbowPlot(Singlet_norm, ndims = 50)
DimHeatmap(Singlet_norm, dims = 1:9, cells = 500, balanced = TRUE, nfeatures = 20)
DimHeatmap(Singlet_norm, dims = 10:18, cells = 500, balanced = TRUE, nfeatures = 20)
DimPlot(Singlet_norm,
        reduction = "pca",
        group.by = "HTO_maxID",
        pt.size = 1,
        cols = my_palette)+
  labs(color = "legend title")+
  ggtitle(label = "PCA")
DimPlot(Singlet_norm, 
        reduction = "pca",
        group.by = "HTO_maxID",
        dims = 3:4,
        pt.size = 1,
        cols = my_palette)+
  labs(color = "legend title")+
  ggtitle(label = "PCA")
### Selcted 11 CP

#UMAP----
Singlet_norm <- RunUMAP(Singlet_norm, dims = 1:11, assay = "SCT")
DimPlot(Singlet_norm, 
        reduction = "umap",
        group.by = "HTO_maxID", 
        pt.size = 1, 
        cols = my_palette)+
  labs(color = "legend title")+
  ggtitle(label = "UMAP 11 CP")

#Clustree-----
Singlet_norm <- FindNeighbors(Singlet_norm, dims = 1:11, verbose = FALSE)
Singlet_norm <- FindClusters(Singlet_norm, 
                             resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 
                                           0.6, 0.7, 0.8, 0.9, 1, 1.1, 1.2))
head(Singlet_norm[[]])
clustree(Singlet_norm, prefix = "SCT_snn_res.")

###   choose resolution 0.3
# Singlet_norm <- FindNeighbors(Singlet_norm, dims = 1:11, verbose = FALSE)
# Singlet_norm <- FindClusters(Singlet_norm, resolution = 0.2)
# Singlet_norm <- RunUMAP(Singlet_norm, dims = 1:11, assay = "SCT")
# 
# DimPlot(Singlet_norm, 
#         reduction = "umap", 
#         pt.size = 1)+
#   ggtitle(label = "UMAP resolution : 0.2")

reso = c(0.3) #, 0.4, 0.5, 0.6, 0.7, 0.8
for (i in reso) {
  Singlet_norm <- FindNeighbors(Singlet_norm, dims = 1:11, verbose = FALSE)
  Singlet_norm <- FindClusters(Singlet_norm, resolution = i)
  Singlet_norm <- RunUMAP(Singlet_norm, dims = 1:11, assay = "SCT")
  
  a<-DimPlot(Singlet_norm, 
          reduction = "umap", 
          pt.size = 1,
          label = T)+
    ggtitle(label = paste0( "UMAP resolution : " , i , " "))
  print(a)
}


#####
condition <- as.data.frame(Singlet_norm@meta.data[["HTO_maxID"]]) #creer un df de tous tes idividus
Singlet_norm[["orig.ident"]] <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  #on passe la 1er colonne
Idents(Singlet_norm) <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  # reviens a faire un levels
levels(Singlet_norm)
####

levels(Singlet_norm)
# test2 <- FindMarkers(Singlet_norm, 
#                      ident.1 = "Ly49p-CD8ab", 
#                      ident.2 = "Ly49p-CD8aa",
#                      min.pct = 0.1, 
#                      logfc.threshold = log(2))
# test3 <- FindMarkers(Singlet_norm, 
#                      ident.1 = "VV-Ly49p-CD8aa", 
#                      ident.2 = "Ly49p-CD8aa",
#                      min.pct = 0.1, 
#                      logfc.threshold = log(2))

Singlet_norm.markers <- FindAllMarkers(Singlet_norm, 
                                       only.pos = F, 
                                       min.pct = 0.1, 
                                       logfc.threshold = log(2))
Singlet_norm.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

Singlet_norm.markers %>%
  group_by(cluster) %>%
  top_n(n = 30, wt = avg_log2FC) -> top10
DoHeatmap(Singlet_norm, features = top10$gene) + NoLegend()

reso = c(0.5)
for (i in reso) {
  Singlet_norm <- FindNeighbors(Singlet_norm, dims = 1:11, verbose = FALSE)
  Singlet_norm <- FindClusters(Singlet_norm, resolution = i)
  Singlet_norm <- RunUMAP(Singlet_norm, dims = 1:11, assay = "SCT")
  
  a<-DimPlot(Singlet_norm, 
             reduction = "umap", 
             pt.size = 1,
             label = T)+
    ggtitle(label = paste0( "UMAP resolution : " , i , " "))
  print(a)
}
levels(Singlet_norm)

Singlet_norm.markers <- FindAllMarkers(Singlet_norm, 
                                       only.pos = F, 
                                       min.pct = 0.1, 
                                       logfc.threshold = log(2))
Singlet_norm.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

Singlet_norm.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> top10
DoHeatmap(Singlet_norm, features = top10$gene) + NoLegend()+
  ggtitle(label = "resolution : 0.5")













