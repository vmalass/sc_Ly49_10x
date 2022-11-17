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
                        keys=inter$gene, 
                        column="ENTREZID", 
                        keytype="SYMBOL", 
                        multiVals="first")

Kres <- enrichKEGG(gene         = geneListEntrez,
                   organism     = "mmu",
                   pvalueCutoff = 0.05)

Kresmat = as.matrix(Kres[])

barplot(Kres, showCategory=20, order=T) + ggtitle("Barplot of functional enrichment by KEGG")
dotplot(Kres, showCategory=20) + ggtitle("Dotplot of functional enrichment by KEGG")


################################################################################
################################################################################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(row.names(DE_5_vs_8), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = DE_5_vs_8[row.names(DE_5_vs_8) %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               nPerm        = 10000,
               minGSSize    = 3,
               maxGSSize    = 800,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2mat = as.matrix(kk2[])
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
ridgeplot(kk2) + labs(x = "enrichment distribution")

################################################################################
################################################################################
################################################################################

# Convert gene IDs for gseKEGG function
# We will lose some genes here because not all IDs will be converted
ids<-bitr(inter$gene, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Mm.eg.db)
# remove duplicate IDS (here I use "ENSEMBL", but it should be whatever was selected as keyType)
dedup_ids = ids[!duplicated(ids[c("SYMBOL")]),]

# Create a new dataframe df2 which has only the genes which were successfully mapped using the bitr function above
df2 = inter[inter$gene %in% dedup_ids$SYMBOL,]

# Create a new column in df2 with the corresponding ENTREZ IDs
df2$Y = dedup_ids$ENTREZID

# Create a vector of the gene unuiverse
kegg_gene_list <- df2$avg_log2FC

# Name vector with ENTREZ ids
names(kegg_gene_list) <- df2$Y

# omit any NA values 
kegg_gene_list<-na.omit(kegg_gene_list)

# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list = sort(kegg_gene_list, decreasing = TRUE)

kegg_organism = "mmu"
kk2 <- gseKEGG(geneList     = kegg_gene_list,
               organism     = kegg_organism,
               minGSSize    = 3,
               pvalueCutoff = 0.05,
               pAdjustMethod = "none",
               keyType       = "ncbi-geneid")

kk2mat = as.matrix(kk2[])
dotplot(kk2, showCategory = 10, title = "Enriched Pathways" , split=".sign") + facet_grid(.~.sign)
dotplot(kk2, showCategory=20) + ggtitle("Dotplot of functional enrichment by GO")

ridgeplot(kk2) + labs(x = "enrichment distribution")

# GO ---------------------------------------------------------------------------
geneListEnsembl= mapIds(org.Mm.eg.db, 
                        keys=inter$gene, 
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
                           keys=inter$gene, 
                           column="ENTREZID", 
                           keytype="SYMBOL", 
                           multiVals="first")

enrichedRes_GS <- enrichPathway(gene = geneListEntrez_PA,
                                organism = 'mouse',
                                pvalueCutoff = 0.05)

enrichedRes_GS_mat = as.matrix(enrichedRes_GS[])

barplot(enrichedRes_GS, showCategory=10, orderby="x") + ggtitle("Barplot of functional enrichment by REACTOM_PA")
dotplot(enrichedRes_GS, showCategory=20) + ggtitle("Dotplot of functional enrichment by REACTOM_PA")















