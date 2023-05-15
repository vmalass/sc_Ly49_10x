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
plan("multiprocess", workers = 9) # activate parallelization

Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_3.RDS")

my_palette = c("Naive"= "#AA4499",
               "Ly49n"= "#332288",
               "VV-B8R"= "#DDCC77",
               "Ly49p-CD8ab"= "#117733",
               "Ly49p-CD8aa"= "#88CCEE",
               "VV-Ly49p-CD8ab"= "#999933",
               "VV-Ly49p-CD8aa"= "#882255",
               "Ly49_Ag_spe" = "red",
               "G1" = "salmon2",
               "G2M" = "green4",
               "S" = "dodgerblue")

# 3 Changement des levels pop---------------------------------------------------
levels(Singlet_norm)
condition <- as.data.frame(Singlet_norm@meta.data[["HTO_maxID"]]) #creer un df de tous tes idividus
Singlet_norm[["orig.ident"]] <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  #on passe la 1er colonne
Idents(Singlet_norm) <- condition$`Singlet_norm@meta.data[["HTO_maxID"]]`  # reviens a faire un levels
levels(Singlet_norm)

# 4 Nouveau objet seurat--------------------------------------------------------
# data_filtre <- subset(Singlet_norm, HTO_maxID != "VV-B8R")
data_filtre <- subset(Singlet_norm, HTO_maxID != "Ly49n")
data_filtre <- subset(data_filtre, HTO_maxID != "Naive")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "3")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "7")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "0")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "2")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "1")
data_filtre <- subset(data_filtre, SCT_snn_res.0.4 != "6")

levels(data_filtre)

### UMAP cluster ###
umap_pop <- DimPlot(data_filtre,
                    reduction = "umap",
                    group.by = "HTO_maxID",
                    pt.size = 1,
                    label = F)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP avec VV-B8R et les Ly49 Ag spé"))

umap_clus <- DimPlot(data_filtre,
                     reduction = "umap",
                     group.by = "SCT_snn_res.0.4",
                     pt.size = 1,
                     label = T)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP avec VV-B8R et les Ly49 Ag spé"))

plot_grid(umap_pop, umap_clus, labels=c("A", "B"), ncol = 1, nrow = 2)  # Plot avec 2 figures

# 5 DE--------------------------------------------------------------------------
## 5.1 VV-B8R vs Ly49 Ag spé----------------------------------------------------
Ly49_ag_spe_VVB8R <- FindMarkers(data_filtre,
                           group.by = "HTO_maxID",
                           ident.1 = c("VV-B8R"),
                           ident.2 = c("Ly49p-CD8ab", "VV-Ly49p-CD8aa", "Ly49p-CD8aa", "VV-Ly49p-CD8ab"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

Ly49_ag_spe_VVB8R <- arrange(Ly49_ag_spe_VVB8R, avg_log2FC)

DoHeatmap(data_filtre, features = rownames(Ly49_ag_spe_VVB8R), group.by = "HTO_maxID") +
  ggtitle(label = paste0("Heatmaps du DE VV-B8R et Ly49 Ag spe"))
DoHeatmap(data_filtre, features = rownames(Ly49_ag_spe_VVB8R), group.by = "SCT_snn_res.0.4") +
  ggtitle(label = paste0("Heatmaps du DE VV-B8R et Ly49 Ag spe"))

# 6 Enrichissement--------------------------------------------------------------
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

## 6.1 VV-B8R vs Ly49 Ag spé----------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(Ly49_ag_spe_VVB8R),
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
                        keys= rownames(Ly49_ag_spe_VVB8R),
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
                           keys= rownames(Ly49_ag_spe_VVB8R), 
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








