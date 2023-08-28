DE entre les differentes pop

# 1 Library---------------------------------------------------------------------
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

if (!require('ggplot2')) install.packages('ggplot2'); library('ggplot2')
if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')
if (!require('cowplot')) install.packages('cowplot'); library('cowplot')
if (!require('org.Mm.eg.db')) BiocManager::install('org.Mm.eg.db'); library('org.Mm.eg.db')  # mouse
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')
if (!require('enrichplot')) BiocManager::install('enrichplot'); library('enrichplot')
if (!require('ReactomePA')) BiocManager::install('ReactomePA'); library('ReactomePA')

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

# 3 DE -------------------------------------------------------------------------
## 3.1 Ly49 aa vs VVaa----------------------------------------------------------

unique(Singlet_norm@meta.data[["orig.ident"]])

DE_aa_VVaa <- FindMarkers(Singlet_norm,
                           group.by = "pop_Ag_spe",
                           ident.1 = c("Ly49p_CD8aa"),
                           ident.2 = c("VV_Ly49p_CD8aa"),
                           only.pos = F,
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                           min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

DE_aa_VVaa <- arrange(DE_aa_VVaa, avg_log2FC)
a <- data.frame(rownames(DE_aa_VVaa))

DoHeatmap(Singlet_norm, features = rownames(DE_aa_VVaa), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8aa vs VV_Ly49p_CD8aa"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.2 Ly49 ab vs VVab----------------------------------------------------------

DE_ab_VVab <- FindMarkers(Singlet_norm,
                          group.by = "pop_Ag_spe",
                          ident.1 = c("Ly49p_CD8ab"),
                          ident.2 = c("VV_Ly49p_CD8ab"),
                          only.pos = F,
                          min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                          logfc.threshold = 0.5,  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49
                          min.diff.pct = 0.25)  # Pre-filter features whose detection percentages across the two groups are similar 

DE_ab_VVab <- arrange(DE_ab_VVab, avg_log2FC)
a <- data.frame(rownames(DE_ab_VVab))

DoHeatmap(Singlet_norm, features = rownames(DE_ab_VVab), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8ab vs VV_Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.3 VV Ly49 aa vs ab----------------------------------------------------------

DE_VV_aa_ab <- FindMarkers(Singlet_norm,
                          group.by = "pop_Ag_spe",
                          ident.1 = c("VV_Ly49p_CD8aa"),
                          ident.2 = c("VV_Ly49p_CD8ab"),
                          only.pos = F,
                          min.diff.pct = 0.25,  # Pre-filter features whose detection percentages across the two groups are similar
                          min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                          logfc.threshold = 0.5)  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49

DE_VV_aa_ab <- arrange(DE_VV_aa_ab, avg_log2FC)
a <- data.frame(rownames(DE_VV_aa_ab))

DoHeatmap(Singlet_norm, features = rownames(DE_VV_aa_ab), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE VV_Ly49p_CD8ab vs VV_Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.4  Ly49 aa vs ab----------------------------------------------------------

DE_aa_ab <- FindMarkers(Singlet_norm,
                           group.by = "pop_Ag_spe",
                           ident.1 = c("Ly49p_CD8aa"),
                           ident.2 = c("Ly49p_CD8ab"),
                           only.pos = F,
                           min.diff.pct = 0.25,  # Pre-filter features whose detection percentages across the two groups are similar
                           min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                           logfc.threshold = 0.5)  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49

DE_aa_ab <- arrange(DE_aa_ab, avg_log2FC)
a <- data.frame(rownames(DE_aa_ab))

DoHeatmap(Singlet_norm, features = rownames(DE_aa_ab), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8aa vs Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.5  Ly49 VVaa&ab vs non VVaa&ab----------------------------------------------------------

DE_aaab_VVaaab <- FindMarkers(Singlet_norm,
                        group.by = "pop_Ag_spe",
                        ident.1 = c("Ly49p_CD8aa", "Ly49p_CD8ab"),
                        ident.2 = c("VV_Ly49p_CD8aa", "VV_Ly49p_CD8ab"),
                        only.pos = F,
                        min.diff.pct = 0.25,  # Pre-filter features whose detection percentages across the two groups are similar
                        min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                        logfc.threshold = 0.5)  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49

DE_aaab_VVaaab <- arrange(DE_aaab_VVaaab, avg_log2FC)
a <- data.frame(rownames(DE_aaab_VVaaab))

DoHeatmap(Singlet_norm, features = rownames(DE_aaab_VVaaab), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8aa vs Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


## 3.4  Ly49 aa vs ab----------------------------------------------------------

DE_aa_ab <- FindMarkers(Singlet_norm,
                        group.by = "pop_Ag_spe",
                        ident.1 = c("Ly49p_CD8aa"),
                        ident.2 = c("Ly49p_CD8ab"),
                        only.pos = F,
                        min.diff.pct = 0.25,  # Pre-filter features whose detection percentages across the two groups are similar
                        min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                        logfc.threshold = 0.5)  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49

DE_aa_ab <- arrange(DE_aa_ab, avg_log2FC)
a <- data.frame(rownames(DE_aa_ab))

DoHeatmap(Singlet_norm, features = rownames(DE_aa_ab), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8aa vs Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 14),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))

## 3.6  Ly49 VVaa&ab vs non VVaa&ab----------------------------------------------------------

DE_ly49_Agspe <- FindMarkers(Singlet_norm,
                              group.by = "pop_Ag_spe",
                              ident.1 = c("Ly49p_CD8aa", "Ly49p_CD8ab", "VV_Ly49p_CD8aa", "VV_Ly49p_CD8ab"),
                              ident.2 = c("Ly49_Ag_spe"),
                              only.pos = F,
                              min.diff.pct = 0.25,  # Pre-filter features whose detection percentages across the two groups are similar
                              min.pct = 0.1,  # Pre-filter features that are detected at <10% frequency in either Ly49 Ag spe ("8", "5", "4") or Ly49 ("0", "1", "2", "6")
                              logfc.threshold = 0.5)  # Pre-filter features that have less than a 1/2-fold change between the average expression of  Ly49 Ag spe vs Ly49

DE_ly49_Agspe <- arrange(DE_ly49_Agspe, avg_log2FC)
a <- data.frame(rownames(DE_ly49_Agspe))

DoHeatmap(Singlet_norm, features = rownames(DE_ly49_Agspe), group.by = "pop_Ag_spe") +
  # NoLegend() +
  ggtitle(label = paste0("Heatmaps du DE Ly49p_CD8aa vs Ly49p_CD8ab"))+ 
  theme(plot.title = element_text(size=14),
        axis.title = element_text(size = 14),
        axis.text = element_text(size = 8),
        legend.title = element_text(size=15),
        legend.text = element_text(size=15))


# 4 Enrichissement--------------------------------------------------------------
### univers Kegg reactom ####
# geneUniverse <- as.data.frame(Singlet_norm@assays[["SCT"]]@counts@Dimnames[[1]])
geneUniverse <- as.data.frame(Singlet_norm@assays[["SCT"]]@var.features)
names(geneUniverse) <- "gene"

geneUniverse = mapIds(org.Mm.eg.db, 
                      keys= geneUniverse$gene,
                      column="ENTREZID", 
                      keytype="SYMBOL", 
                      multiVals="first")
geneUniverse <- na.omit(geneUniverse)  # pas de mappage de certain gene
print(paste0("Pourcentage de gène identifié Kegg et Reactom (NA) : ", round(100 - (length(geneUniverse)/length(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])*100), 2), " %"))

### univers GO ###
# geneUniverse_GO <- as.data.frame(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])
geneUniverse_GO <- as.data.frame(Singlet_norm@assays[["SCT"]]@var.features)
names(geneUniverse_GO) <- "gene"

geneUniverse_GO <- mapIds(org.Mm.eg.db, 
                          keys= geneUniverse_GO$gene,
                          column="ENSEMBL", 
                          keytype="SYMBOL", 
                          multiVals="first")
geneUniverse_GO <- na.omit(geneUniverse_GO)  # pas de mappage de certain gene
print(paste0("Pourcentage de gène identifié GO (NA) : ", round(100 - (length(geneUniverse_GO)/length(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])*100), 2), " %"))

## 6.1 VV-B8R vs Ly49 Ag spé----------------------------------------------------
# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys= rownames(DE_aaab_VVaaab),
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
                        keys= rownames(DE_aaab_VVaaab),
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
                           keys= rownames(DE_aaab_VVaaab), 
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

plot_grid(go, pa, labels=c( "Go", "Reactom_PA"), ncol = 2, nrow = 1)  # Plot avec 3 figures

