

# 1 Library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')








# 2 import data-----------------------------------------------------------------
rm(list = ls())
plan("multiprocess", workers = 9) # activate parallelization
Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_5.RDS")

my_palette = c("Naive"= "#AA4499",
               "Ly49n"= "#332288",
               "VV-B8R"= "#DDCC77",
               "VV_B8R"= "#DDCC77",
               "Ly49p-CD8ab"= "#117733",
               "Ly49p_CD8ab"= "#117733",
               "Ly49p-CD8aa"= "#88CCEE",
               "Ly49p_CD8aa"= "#88CCEE",
               "VV-Ly49p-CD8ab"= "#999933",
               "VV_Ly49p_CD8ab"= "#999933",
               "VV-Ly49p-CD8aa"= "#882255",
               "VV_Ly49p_CD8aa"= "#882255",
               "Ly49_Ag_spe" = "red",
               "G1" = "salmon2",
               "G2M" = "green4",
               "S" = "dodgerblue")

# 3 PCA UMAP--------------------------------------------------------------------
levels(Singlet_norm)
### PCA ####
pca_pop <- DimPlot(Singlet_norm,
                   reduction = "pca",
                   group.by = "pop_Ag_spe",
                   pt.size = 1,
                   cols = my_palette,
                   label = T)+
  labs(color = "legend :")+
  ggtitle(label = "PCA population")

pca_clus <- DimPlot(Singlet_norm,
                    reduction = "pca",
                    group.by = "SCT_snn_res.0.4",
                    pt.size = 1,
                    label = T)+
  labs(color = "legend :")+
  ggtitle(label = "PCA cluster")

### UMAP ###
umap_pop <- DimPlot(Singlet_norm,
                    reduction = "umap",
                    group.by = "pop_Ag_spe",
                    pt.size = 1,
                    cols = my_palette,
                    label = T) +
  labs(color = "Legend :") +
  ggtitle(label = paste0("UMAP avec les 11 premieres CP"))

### UMAP cluster ###
umap_clus <- DimPlot(Singlet_norm,
                     reduction = "umap",
                     group.by = "SCT_snn_res.0.4",
                     pt.size = 1,
                     label = T)+
  labs(color = "legend ")+
  ggtitle(label = paste0( "UMAP avec les 11 premieres 11 CP et resolution : 0.4"))


plot_grid(pca_pop, pca_clus, umap_pop, umap_clus, labels=c("A", "B", "C", "D"), ncol = 2, nrow = 2)  # Plot avec 2 figures

# 4 Changement des levels cluster-----------------------------------------------
levels(Singlet_norm)
condition <- as.data.frame(Singlet_norm@meta.data[["pop_Ag_spe"]]) #creer un df de tous tes idividus
Singlet_norm[["orig.ident"]] <- condition$`Singlet_norm@meta.data[["pop_Ag_spe"]]`  #on passe la 1er colonne
Idents(Singlet_norm) <- condition$`Singlet_norm@meta.data[["pop_Ag_spe"]]`  # reviens a faire un levels
levels(Singlet_norm)

# 5 DE--------------------------------------------------------------------------
# Ly49_Ag_spe vs Ly49 non VV
Ly49_ag_spe <- FindMarkers(Singlet_norm,
                           # group.by = "pop_Ag_spe",
                           ident.1 = c("Ly49_Ag_spe"),
                           ident.2 = c("Ly49p_CD8ab", "Ly49p_CD8aa", "VV_Ly49p_CD8ab", "VV_Ly49p_CD8aa"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_ag_spe <- arrange(Ly49_ag_spe, avg_log2FC)
a <- data.frame(rownames(Ly49_ag_spe))

DoHeatmap(Singlet_norm, features = rownames(Ly49_ag_spe), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49_Ag_spe et Ly49p_CD8ab Ly49p_CD8aa"))+ 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 5),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))



# Ly49_Ag_spe vs Ly49 non VV
Ly49_ag_spe <- FindMarkers(Singlet_norm,
                           # group.by = "pop_Ag_spe",
                           ident.1 = c("Ly49_Ag_spe"),
                           ident.2 = c("VV_B8R"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_ag_spe <- arrange(Ly49_ag_spe, avg_log2FC)
a <- data.frame(rownames(Ly49_ag_spe))

data_filtre <- subset(Singlet_norm, pop_Ag_spe == c("Ly49_Ag_spe", "VV_B8R"))

DoHeatmap(data_filtre, features = rownames(Ly49_ag_spe), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49_Ag_spe vs VV-B8R"))+ 
  theme(plot.title = element_text(size=24),
        axis.title = element_text(size = 24),
        axis.text = element_text(size = 12),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


write.csv(a, "list_DE_Ly49AgSpe_VVB8R.csv")

