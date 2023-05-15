# 1 Library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('org.Mm.eg.db')) BiocManager::install('org.Mm.eg.db'); library('org.Mm.eg.db')  # mouse
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')
if (!require('enrichplot')) BiocManager::install('enrichplot'); library('enrichplot')
if (!require('ReactomePA')) BiocManager::install('ReactomePA'); library('ReactomePA')
if (!require('stringr')) install.packages('stringr'); library('stringr')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')


library(cowplot)


# 2 import data-----------------------------------------------------------------
rm(list = ls())
# plan("multiprocess", workers = 9) # activate parallelization

Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_3.RDS")

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
### PCA ####
DimPlot(Singlet_norm,
        reduction = "pca",
        group.by = "HTO_maxID",
        pt.size = 1,
        cols = my_palette)+
  labs(color = "legend ")+
  ggtitle(label = "PCA")

### UMAP ###
DimPlot(Singlet_norm,
        reduction = "umap",
        group.by = "HTO_maxID",
        pt.size = 1,
        cols = my_palette,
        label = T) +
  labs(color = "Legend :") +
  ggtitle(label = paste0("UMAP avec les 11 premieres CP"))

### UMAP cluster ###
DimPlot(Singlet_norm,
        reduction = "umap",
        group.by = "SCT_snn_res.0.4",
        pt.size = 1,
        label = T)+
  labs(color = "legend ")+
  ggtitle(label = paste0( "UMAP avec les 11 premieres 11 CP et resolution : 0.4"))

# 4 Changement des levels pop---------------------------------------------------
levels(Singlet_norm)
condition <- as.data.frame(Singlet_norm@meta.data[["HTO_maxID"]]) #creer un df de tous tes idividus
Singlet_norm[["orig.ident"]] <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  #on passe la 1er colonne
Idents(Singlet_norm) <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  # reviens a faire un levels
levels(Singlet_norm)

# 5 Nouveau objet seurat--------------------------------------------------------
data_filtre <- subset(Singlet_norm, HTO_maxID != "VV-B8R")
data_filtre <- subset(data_filtre, HTO_maxID != "Ly49n")
data_filtre <- subset(data_filtre, HTO_maxID != "Naive")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "3")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "7")

levels(data_filtre)
### UMAP cluster ###
umap_pop <- DimPlot(data_filtre,
        reduction = "umap",
        group.by = "HTO_maxID",
        pt.size = 1,
        label = F)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP sans les VV-B8R, Ly49n, Naive et cluster 3, 7"))

umap_clus <- DimPlot(data_filtre,
        reduction = "umap",
        group.by = "SCT_snn_res.0.4",
        pt.size = 1,
        label = T)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP sans les VV-B8R, Ly49n, Naive et cluster 3, 7"))

plot_grid(umap_pop, umap_clus, labels=c("A", "B"), ncol = 2, nrow = 1)  # Plot avec 2 figures



# 6 Changement des levels cluster-----------------------------------------------
levels(data_filtre)
condition <- as.data.frame(data_filtre@meta.data[["SCT_snn_res.0.4"]]) #creer un df de tous tes idividus
data_filtre[["orig.ident"]] <- condition$`data_filtre@meta.data[["SCT_snn_res.0.4"]]`  #on passe la 1er colonne
Idents(data_filtre) <- condition$`data_filtre@meta.data[["SCT_snn_res.0.4"]]`  # reviens a faire un levels
levels(data_filtre)

# 7 identification des pop dans les cluster VV-B8R------------------------------
f <- subset(data_filtre, idents = c("8", "5", "4"))
metadata_f <- as.data.frame(f@meta.data[["HTO_maxID"]])
metadata_f <- cbind(metadata_f, as.data.frame(f@meta.data[["SCT_snn_res.0.4"]]))
colnames(metadata_f)<- c("population","cluster")

metadata_f %>%
  count(population, cluster) %>%
  arrange(cluster, population) %>%
  mutate(pct_total = n / sum(n) * 100) -> e

print(e)

# 8 DE--------------------------------------------------------------------------
## 8.1 Ly49 Ag spé vs Ly49------------------------------------------------------
Ly49_ag_spe <- FindMarkers(data_filtre,
                           group.by = "SCT_snn_res.0.4",
                           ident.1 = c("8", "5", "4"),
                           ident.2 = c("0", "1", "2", "6"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_ag_spe <- arrange(Ly49_ag_spe, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_ag_spe), group.by = "SCT_snn_res.0.4") +
  NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49 Ag spe (8, 5, 4) et Ly49 (0, 1, 2 , 6)"))


## 8.2 Ly49 VV vs Ly49----------------------------------------------------------
Ly49_VV <- FindMarkers(data_filtre,
                           group.by = "SCT_snn_res.0.4",
                           ident.1 = c("0", "2"),
                           ident.2 = c("1"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_VV <- arrange(Ly49_VV, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_VV), group.by = "SCT_snn_res.0.4") +
  NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49 VV (0, 2) et Ly49 (1)"))


## 8.3 Ly49 IFN vs Ly49----------------------------------------------------------
Ly49_IFN <- FindMarkers(data_filtre,
                       group.by = "SCT_snn_res.0.4",
                       ident.1 = c("6"),
                       ident.2 = c("1"),
                       only.pos = F,
                       min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                       logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                       min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_IFN <- arrange(Ly49_IFN, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_IFN), group.by = "SCT_snn_res.0.4") +
  NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49 IFN (6) et Ly49 (1)"))


## 8.4 Ly49_CD8aa VV vs Ly49 IFN------------------------------------------------
Ly49_CD8aa_VV <- FindMarkers(data_filtre,
                        group.by = "SCT_snn_res.0.4",
                        ident.1 = c("2"),
                        ident.2 = c("6"),
                        only.pos = F,
                        min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                        logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                        min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_CD8aa_VV <- arrange(Ly49_CD8aa_VV, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_IFN), group.by = "SCT_snn_res.0.4") +
  NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49_CD8aa VV (2) et Ly49 IFN (6)"))

## 8.5 Ly49 Ag spé vs Ly49 non VV-----------------------------------------------
Ly49_ag_spe_Ly48_nonVV <- FindMarkers(data_filtre,
                           group.by = "SCT_snn_res.0.4",
                           ident.1 = c("8", "5", "4"),
                           ident.2 = c("1"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_ag_spe_Ly48_nonVV <- arrange(Ly49_ag_spe_Ly48_nonVV, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_ag_spe_Ly48_nonVV), group.by = "SCT_snn_res.0.4") +
  NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49 Ag spe (8, 5, 4) et Ly49 non VV (1)"))+
  theme(axis.text.y = element_text(size = 6))


# 9 Enrichissement--------------------------------------------------------------
### univers Kegg reactom ####
geneUniverse <- as.data.frame(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])
geneUniverse = mapIds(org.Mm.eg.db, 
                        keys= geneUniverse$`Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]]`,
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")
geneUniverse <- na.omit(geneUniverse)  # pas de mappage de certain gene
print(paste0("Pourcentage de gène non identifié Kegg et Reactom (NA) : ", round(100 - (length(geneUniverse)/length(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])*100), 2), " %"))

### univers GO ###
geneUniverse_GO <- as.data.frame(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])
geneUniverse_GO <- mapIds(org.Mm.eg.db, 
                        keys= geneUniverse_GO$`Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]]`,
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")
geneUniverse_GO <- na.omit(geneUniverse_GO)  # pas de mappage de certain gene
print(paste0("Pourcentage de gène non identifié GO (NA) : ", round(100 - (length(geneUniverse_GO)/length(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])*100), 2), " %"))


## 9.1 Ly49 VV vs Ly49----------------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_VV),
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "mmu",
                   universe = geneUniverse,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
kegg <- dotplot(Kres, showCategory=10) + ggtitle("Dotplot of functional enrichment by KEGG")

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_VV),
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  universe = geneUniverse_GO,
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by GO")
go <- dotplot(enrichedRes_GO, showCategory=10) + ggtitle("Dotplot of functional enrichment by GO")

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys= rownames(Ly49_VV), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_PA <- enrichPathway(gene          = geneListEntrez_PA,
                                universe = geneUniverse,
                                organism      = 'mouse',
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

enrichedRes_PA_mat = as.matrix(enrichedRes_PA[])

barplot(enrichedRes_PA, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
pa <- dotplot(enrichedRes_PA, showCategory=10) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")


plot_grid(kegg, go, pa, labels=c("Kegg", "Go", "Reactom_PA"), ncol = 3, nrow = 1)  # Plot avec 3 figures

## 9.2 Ly49 Ag spe vs Ly49------------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_ag_spe),
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "mmu",
                   universe = geneUniverse,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
kegg <- dotplot(Kres, showCategory=10) + ggtitle("Dotplot of functional enrichment by KEGG")

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_ag_spe),
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  universe = geneUniverse_GO,
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by GO")
go <- dotplot(enrichedRes_GO, showCategory=10) + ggtitle("Dotplot of functional enrichment by GO")

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys= rownames(Ly49_ag_spe), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_PA <- enrichPathway(gene          = geneListEntrez_PA,
                                organism      = 'mouse',
                                universe = geneUniverse,
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

enrichedRes_PA_mat = as.matrix(enrichedRes_PA[])

barplot(enrichedRes_PA, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
pa <- dotplot(enrichedRes_PA, showCategory=10) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")


plot_grid(kegg, go, pa, labels=c("Kegg", "Go", "Reactom_PA"), ncol = 3, nrow = 1)  # Plot avec 3 figures

## 9.3 Ly49 IFN vs Ly49------------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_IFN),
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "mmu",
                   universe = geneUniverse,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
kegg <- dotplot(Kres, showCategory=10) 

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_IFN),
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  universe = geneUniverse_GO,
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by GO")
go <- dotplot(enrichedRes_GO, showCategory=10) 

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys= rownames(Ly49_IFN), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_PA <- enrichPathway(gene          = geneListEntrez_PA,
                                organism      = 'mouse',
                                universe = geneUniverse,
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

enrichedRes_PA_mat = as.matrix(enrichedRes_PA[])

barplot(enrichedRes_PA, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
pa <- dotplot(enrichedRes_PA, showCategory=10) 

plot_grid( go, pa, labels=c("Go", "Reactom_PA"), ncol = 3, nrow = 1)  # Plot avec 3 figures

## 9.4 Ly49 CD8aa VV vs Ly49 IFN------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_CD8aa_VV),
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "mmu",
                   universe = geneUniverse,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
kegg <- dotplot(Kres, showCategory=10) 

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_CD8aa_VV),
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  universe = geneUniverse_GO,
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by GO")
go <- dotplot(enrichedRes_GO, showCategory=10) 

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys= rownames(Ly49_CD8aa_VV), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_PA <- enrichPathway(gene          = geneListEntrez_PA,
                                organism      = 'mouse',
                                universe = geneUniverse,
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

enrichedRes_PA_mat = as.matrix(enrichedRes_PA[])

barplot(enrichedRes_PA, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
pa <- dotplot(enrichedRes_PA, showCategory=10) 

plot_grid(kegg, go, pa, labels=c("Kegg", "Go", "Reactom_PA"), ncol = 3, nrow = 1)  # Plot avec 3 figures


## 9.5 Ly49 Ag spe vs Ly49 non VV------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_ag_spe_Ly48_nonVV),
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene          = geneListEntrez,
                   organism      = "mmu",
                   universe = geneUniverse,
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05, # slider
                   qvalueCutoff  = 0.05) # slider

Kresmat = as.matrix(Kres[])

# barplot(Kres, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
kegg <- dotplot(Kres, showCategory=10) 

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_ag_spe_Ly48_nonVV),
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  universe = geneUniverse_GO,
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05, # slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

# barplot(enrichedRes_GO, showCategory=10, order=T) + ggtitle("Barplot of functional enrichment by GO")
go <- dotplot(enrichedRes_GO, showCategory=10) 

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys= rownames(Ly49_ag_spe_Ly48_nonVV), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_PA <- enrichPathway(gene          = geneListEntrez_PA,
                                organism      = 'mouse',
                                universe = geneUniverse,
                                pAdjustMethod = "BH", 
                                pvalueCutoff  = 0.05, # slider
                                qvalueCutoff  = 0.05) # slider

enrichedRes_PA_mat = as.matrix(enrichedRes_PA[])

# barplot(enrichedRes_PA, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
pa <- dotplot(enrichedRes_PA, showCategory=10) 

plot_grid(kegg, go, pa, labels=c("Kegg", "Go", "Reactom_PA"), ncol = 3, nrow = 1)  # Plot avec 3 figures








a <- setdiff(rownames(Ly49_ag_spe_Ly48_nonVV), rownames(Ly49_VV))
b <- setdiff(rownames(Ly49_VV), rownames(Ly49_ag_spe_Ly48_nonVV))
c <- intersect(rownames(Ly49_VV), rownames(Ly49_ag_spe_Ly48_nonVV))
d <- union(rownames(Ly49_VV), rownames(Ly49_ag_spe_Ly48_nonVV))












