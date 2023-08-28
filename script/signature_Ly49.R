if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
if (!require('openxlsx')) install.packages('openxlsx'); library('openxlsx')
if (!require('dplyr')) install.packages('dplyr'); library('dplyr')

# Data import-------------------------------------------------------------------
Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_5.RDS")

# Visualization-----------------------------------------------------------------
DimPlot(Singlet_norm,
        reduction = "umap",
        group.by = "SCT_snn_res.0.4", 
        label = T)

# DE for Ly49 signature---------------------------------------------------------
DE_naive_Ly49 <- FindMarkers(Singlet_norm,
            group.by = "SCT_snn_res.0.4",
            ident.1 = 3,
            ident.2 = c(1,0,2,6),
            only.pos = F,
            min.pct = 0.1, #  diff de 0.1 entre les deux conditions
            logfc.threshold = 0.5)

DE_naive_Ly49naive <- FindMarkers(Singlet_norm,
                                  group.by = "SCT_snn_res.0.4",
                                  ident.1 = 3,
                                  ident.2 = c(1),
                                  only.pos = F,
                                  min.pct = 0.1, #  diff de 0.1 entre les deux conditions
                                  logfc.threshold = 0.5)

DE_naive_Ly49VV <- FindMarkers(Singlet_norm,
                                  group.by = "SCT_snn_res.0.4",
                                  ident.1 = 3,
                                  ident.2 = c(0,2),
                                  only.pos = F,
                                  min.pct = 0.1, #  diff de 0.1 entre les deux conditions
                                  logfc.threshold = 0.5)


length(intersect(rownames(DE_naive_Ly49), rownames(DE_naive_Ly49naive)))
length(intersect(rownames(DE_naive_Ly49), rownames(DE_naive_Ly49VV)))
length(intersect(rownames(DE_naive_Ly49naive), rownames(DE_naive_Ly49VV)))

## Add column of gene names-----------------------------------------------------
DE_naive_Ly49 %>% 
  arrange(-avg_log2FC) %>% 
  mutate(rownames(DE_naive_Ly49)) %>% 
  rename(gene = `rownames(DE_naive_Ly49)`)-> DE_naive_Ly49

DE_naive_Ly49naive %>% 
  arrange(-avg_log2FC) %>% 
  mutate(rownames(DE_naive_Ly49naive)) %>% 
  rename(gene = `rownames(DE_naive_Ly49naive)`) -> DE_naive_Ly49naive

DE_naive_Ly49VV %>% 
  arrange(-avg_log2FC) %>% 
  mutate(rownames(DE_naive_Ly49VV)) %>% 
  rename(gene = `rownames(DE_naive_Ly49VV)`) -> DE_naive_Ly49VV

# Save data---------------------------------------------------------------------
write.xlsx(DE_naive_Ly49, "SC_Ly49_10X/results/DE_naive_Ly49.xlsx")
write.xlsx(DE_naive_Ly49naive, "SC_Ly49_10X/results/DE_naive_Ly49naive.xlsx")
write.xlsx(DE_naive_Ly49VV, "SC_Ly49_10X/results/DE_naive_Ly49VV.xlsx")
