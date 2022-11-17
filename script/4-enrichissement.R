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
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require('org.Mm.eg.db')) BiocManager::install('org.Mm.eg.db'); library('org.Mm.eg.db')
if (!require('AnnotationDbi')) BiocManager::install('AnnotationDbi'); library('AnnotationDbi')
if (!require('enrichplot')) BiocManager::install('enrichplot'); library('enrichplot')
if (!require('clusterProfiler')) BiocManager::install('clusterProfiler'); library('clusterProfiler')
if (!require('ReactomePA')) BiocManager::install('ReactomePA'); library('ReactomePA')


# KEGG -------------------------------------------------------------------------
geneListEntrez = mapIds(org.Mm.eg.db, 
                        keys=row.names(DE_0_vs_7), 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene         = geneListEntrez,
                   organism     = "mmu",
                   pvalueCutoff = 0.05)

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
dotplot(Kres, showCategory=20) + ggtitle("Dotplot of functional enrichment by KEGG")

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys=row.names(DE_0_vs_7), 
                        column="ENSEMBL", 
                        keytype="SYMBOL", 
                        multiVals="first")

enrichedRes_GO  <- enrichGO(
  gene          = geneListEnsembl,
  keyType       = "ENSEMBL",
  OrgDb         = org.Mm.eg.db,
  ont           = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.05, # slider
  qvalueCutoff  = 0.05,# slider
  readable      = TRUE)

enrichedRes_GO_mat = as.matrix(enrichedRes_GO[])

barplot(enrichedRes_GO, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by GO")
dotplot(enrichedRes_GO, showCategory=20) + ggtitle("Dotplot of functional enrichment by GO")

# REACTOM_PA--------------------------------------------------------------------
geneListEntrez_PA = mapIds(org.Mm.eg.db, 
                           keys=row.names(DE_0_vs_7), 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_GS <- enrichPathway(gene = geneListEntrez_PA,
                                organism = 'mouse',
                                pvalueCutoff = 0.05)

barplot(enrichedRes_GS, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
dotplot(enrichedRes_GS, showCategory=20) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")















