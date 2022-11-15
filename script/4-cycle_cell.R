# 1 Library-----------------------------------------------------------------------
library(Seurat)
library(dplyr)
library(khroma)
library(future) # parallelization in Seurat (NormalizeData() / ScaleData() / 
                # JackStraw() / FindMarkers() / FindIntegrationAnchors() / 
                # FindClusters() - if clustering over multiple resolutions)
library(clustree)
library(ggplot2)
library(openxlsx)
library(tidyr)

# 2 import data-----------------------------------------------------------------
rm(list = ls())
plan("multiprocess", workers = 9) # activate parallelization

Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_3.RDS")

my_palette = c("Naive"= "#AA4499",
               "Ly49n"= "#332288",
               "VV-B8R"= "#DDCC77",
               "Ly49p-CD8ab"= "#117733",
               "Ly49p-CD8aa"= "#88CCEE",
               "VV-Ly49p-CD8ab"= "#999933",
               "VV-Ly49p-CD8aa"= "#882255",
               "G1" = "salmon2",
               "G2M" = "green4",
               "S" = "dodgerblue")

# 3 Classify cells as G2/M, S, or alternatively G1------------------------------

# A list of cell cycle markers, from Tirosh et al, 2015, is loaded with Seurat.  We can
# segregate this list into markers of G2/M phase and markers of S phase
s.genes <- cc.genes$s.genes
s.genes <- tolower(s.genes)
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
s.genes <- firstup(s.genes)
g2m.genes <- cc.genes$g2m.genes
g2m.genes <- tolower(g2m.genes)
g2m.genes <- firstup(g2m.genes)


Singlet_norm <- CellCycleScoring(Singlet_norm,
                               s.features = s.genes,
                               g2m.features = g2m.genes,
                               set.ident = TRUE)
head(Singlet_norm[[]])

## 3.1 Visualization------------------------------------------------------------
DimPlot(Singlet_norm, 
        reduction = "umap", 
        group.by = "Phase",
        pt.size = 1,
        cols = my_palette) +
  ggtitle(label = "Cell Cycle Phase")

# Visualize the distribution of cell cycle markers across
RidgePlot(Singlet_norm, features = c("Pcna", "Top2a", "Mcm6", "Mki67"), ncol = 2, cols = my_palette)

# 4 Save seurat object----------------------------------------------------------
saveRDS(Singlet_norm, "/Users/victor/Documents/JM/singelcell_LY49/data/Singlet_norm_4.RDS")

# 5 Creat metadata--------------------------------------------------------------
metadata <- as.data.frame(Singlet_norm@meta.data[["orig.ident"]])
metadata <- cbind(metadata,as.data.frame(Singlet_norm@meta.data[["Phase"]]))
metadata <- cbind(metadata, as.data.frame(Singlet_norm@meta.data[["SCT_snn_res.0.5"]]))
colnames(metadata)<- c("population","phase","cluster")

## 5.1 Visualization Cell cycle phase by cluster--------------------------------
metadata %>% 
  count(cluster, phase) %>% 
  arrange(cluster, phase) -> a

ggplot(a, aes(fill=phase, y=n, x=cluster)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values= my_palette) +
  theme_minimal() +
  scale_y_continuous(limits = c(0,1500)) +
  ggtitle(label = "Cell cycle phase by cluster")

# clus <- seq(0,10,1)
# for (i in clus) {
#   inter <- a[which(a$cluster %in% i),] 
#   plt <- ggplot(inter, aes(fill=phase, y=n, x=phase)) +
#     geom_bar(position='dodge', stat='identity') +
#     scale_fill_manual(values= my_palette) +
#     theme_minimal() +
#     scale_y_continuous(limits = c(0,1500)) +
#     ggtitle(label = paste0("Cell cycle phase for cluster ", i))
#   print(plt)
# }

ggplot(a, aes(fill=phase, y=n, x=phase)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values= my_palette) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey")) +
  scale_y_continuous(limits = c(0,1500)) +
  ggtitle(label = paste0("Cell cycle phase by cluster ")) +
  facet_wrap(~cluster, scales='free_x')

## 5.2 Visualization Cell cycle phase by cluster--------------------------------
metadata %>% 
  count(population, cluster) %>% 
  arrange(cluster, population) -> a

ggplot(a, aes(fill=population, y=n, x=cluster)) +
  geom_bar(position='dodge', stat='identity') +
  scale_fill_manual(values= my_palette) +
  theme_minimal() +
  scale_y_continuous(limits = c(0,1200)) +
  ggtitle(label = "Population by cluster")

# clus <- seq(0,10,1)
# for (i in clus) {
#   inter <- a[which(a$cluster %in% i),] 
#   plt <- ggplot(inter, aes(fill=population, y=n, x=population)) +
#     geom_bar(position='dodge', stat='identity') +
#     theme_minimal() +
#     scale_fill_manual(values= my_palette) +
#     scale_y_continuous(limits = c(0,1200)) +
#     ggtitle(label = paste0("Population cluster ",i)) +
#     facet_wrap(~population, scales='free_x')
#   print(plt)
# }

ggplot(a, aes(fill=population, y=n, x=population)) +
  geom_bar(position='dodge', stat='identity') +
  theme(axis.text.x = element_text(angle = 45, vjust = 1.05, hjust=1),
        panel.background = element_blank(),
        axis.line = element_line(colour = "grey")) +
  scale_fill_manual(values= my_palette) +
  scale_y_continuous(limits = c(0,1200)) +
  ggtitle(label = paste0("Population by cluster")) +
  facet_wrap(~cluster, scales='free_x')

