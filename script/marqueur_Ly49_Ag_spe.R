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

# 3 creat new classification----------------------------------------------------
a <- as.data.frame(Singlet_norm@assays[["RNA"]]@counts@Dimnames[[1]])

data <- data.frame(Singlet_norm@meta.data[["HTO_maxID"]], Singlet_norm@meta.data[["SCT_snn_res.0.4"]])
colnames(data) <- c("pop", "cluster")
rownames(data) <- colnames(Singlet_norm) # meme nom entre row df et col obj seurat

Ly49_Ag_spe <- data[data$pop == "Ly49p-CD8aa" & c(data$cluster == "4" | data$cluster == "5" | data$cluster == "8") ,] 
Ly49_Ag_spe <- rbind(Ly49_Ag_spe, data[data$pop == "Ly49p-CD8ab" & c(data$cluster == "4" | data$cluster == "5" | data$cluster == "8") ,])
Ly49_Ag_spe <- rbind(Ly49_Ag_spe, data[data$pop == "VV-Ly49p-CD8aa" & c(data$cluster == "4" | data$cluster == "5" | data$cluster == "8") ,])
Ly49_Ag_spe <- rbind(Ly49_Ag_spe, data[data$pop == "VV-Ly49p-CD8ab" & c(data$cluster == "4" | data$cluster == "5" | data$cluster == "8") ,])
Ly49_Ag_spe <- dplyr::mutate(Ly49_Ag_spe, pop2 = "Ly49_Ag_spe")

VV_B8R<- data[data$pop == "VV-B8R" ,]
VV_B8R <- dplyr::mutate(VV_B8R, pop2 = "VV_B8R")
data_fi <- rbind(Ly49_Ag_spe, VV_B8R)

Naive <- data[data$pop == "Naive" ,]
Naive <- dplyr::mutate(Naive, pop2 = "Naive")
data_fi <- rbind(data_fi, Naive)

Ly49n <- data[data$pop == "Ly49n" ,]
Ly49n <- dplyr::mutate(Ly49n, pop2 = "Ly49n")
data_fi <- rbind(data_fi, Ly49n)

Ly49p_CD8aa <- data[data$pop == "Ly49p-CD8aa" & c(data$cluster == "0" | 
                                                    data$cluster == "1" | 
                                                    data$cluster == "2" |
                                                    data$cluster == "6" |
                                                    data$cluster == "7"),] 
Ly49p_CD8aa <- dplyr::mutate(Ly49p_CD8aa, pop2 = "Ly49p_CD8aa")
data_fi <- rbind(data_fi, Ly49p_CD8aa)

Ly49p_CD8ab <- data[data$pop == "Ly49p-CD8ab" & c(data$cluster == "0" | 
                                                    data$cluster == "1" | 
                                                    data$cluster == "2" |
                                                    data$cluster == "3" |
                                                    data$cluster == "6" |
                                                    data$cluster == "7"),] 
Ly49p_CD8ab <- dplyr::mutate(Ly49p_CD8ab, pop2 = "Ly49p_CD8ab")
data_fi <- rbind(data_fi, Ly49p_CD8ab)

VV_Ly49p_CD8aa <- data[data$pop == "VV-Ly49p-CD8aa" & c(data$cluster == "0" | 
                                                          data$cluster == "1" | 
                                                          data$cluster == "2" |
                                                          data$cluster == "6" |
                                                          data$cluster == "7"),] 
VV_Ly49p_CD8aa <- dplyr::mutate(VV_Ly49p_CD8aa, pop2 = "VV_Ly49p_CD8aa")
data_fi <- rbind(data_fi, VV_Ly49p_CD8aa)



VV_Ly49p_CD8ab <- data[data$pop == "VV-Ly49p-CD8ab" & c(data$cluster == "0" | 
                                                          data$cluster == "1" | 
                                                          data$cluster == "2" |
                                                          data$cluster == "3" |
                                                          data$cluster == "6" |
                                                          data$cluster == "7"),] 
VV_Ly49p_CD8ab <- dplyr::mutate(VV_Ly49p_CD8ab, pop2 = "VV_Ly49p_CD8ab")
data_fi <- rbind(data_fi, VV_Ly49p_CD8ab)

### vérif bien les meme nom ###
all(rownames(data_fi) %in% rownames(data))
all(rownames(data_fi) == rownames(data))

### tri pour avoir le même ordre que obj seurat ###
data_fi <- data_fi[colnames(Singlet_norm),]
all(rownames(data_fi) == rownames(data))

Singlet_norm@meta.data$pop_Ag_spe <- data_fi$pop2  # new groupe dans obj seurat 

b <- Singlet_norm@assays$SCT

### Save ###
# saveRDS(Singlet_norm, "~/Documents/JM/singelcell_LY49/data/Singlet_norm_5.RDS")


# 4 Visualisation des marqueurs-------------------------------------------------
umap_pop <- DimPlot(Singlet_norm,
        reduction = "umap",
        group.by = "pop_Ag_spe",
        pt.size = 1,
        label = F,
        cols = my_palette)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP avec les populations"))

umap_clus <- DimPlot(Singlet_norm,
                    reduction = "umap",
                    group.by = "SCT_snn_res.0.4",
                    pt.size = 1,
                    label = T)+
  labs(color = "legend ")+
  ggtitle(label = paste0("UMAP 11 CP avec les cluster"))

plot_grid(umap_pop, umap_clus, labels=c("A", "B"), ncol = 2, nrow = 1)  # Plot avec 2 figures


var <- c("Gzma", "Gzmb", "Gzmk", "Gzmm")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Fcer1g", "Fcgr2b", "Fcgr3", "Fcgrt", "Fcrl1")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Il18r1", "Il18rap", "Il2rb", "Il4ra", "Il6ra", "Il6st", "Il7r")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Cx3cr1", "Cxcr5", "Cxcr6")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Ccl3", "Ccl4", "Ccl5", "Xcl1")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Bcl2", "Bcl2a1b", "Bcl2l11")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Btg1", "Btg2")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Tnf", "Tnfaip3", "Tnfrsf23", "Tnfrsf25", "Tnfrsf26", "Tnfrsf9", 
         "Tnfsf11","Tnfsf4", "Tnfsf8" )
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()


var <- c("Gzma", "Gzmb", "Gzmk", "Gzmm", "Fcer1g", "Fcgr2b", "Fcgr3", "Fcgrt", 
         "Fcrl1", "Il18r1", "Il18rap", "Il2rb")
var2 <- c("Il4ra", "Il6ra", "Il6st", "Il7r", "Cx3cr1", "Cxcr5", "Cxcr6", "Ccl3",
          "Ccl4", "Ccl5", "Xcl1", "Bcl2")
var3 <- c("Bcl2a1b", "Bcl2l11", "Btg1", "Btg2")

VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")
VlnPlot(Singlet_norm, features = var2, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")
VlnPlot(Singlet_norm, features = var3, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")



var <- c("Gzma", "Gzmb", "Gzmk", "Gzmm", "Fcer1g", "Fcgr2b", "Fcgr3", "Fcgrt", 
         "Fcrl1", "Il18r1", "Il18rap", "Il2rb", "Il4ra", "Il6ra", "Il6st", 
         "Il7r", "Cx3cr1", "Cxcr5", "Cxcr6", "Ccl3", "Ccl4", "Ccl5", "Xcl1",
         "Bcl2", "Bcl2a1b", "Bcl2l11", "Btg1", "Btg2", "Tnf", "Tnfaip3", 
         "Tnfrsf23", "Tnfrsf25", "Tnfrsf26", "Tnfrsf9", "Tnfsf11","Tnfsf4", 
         "Tnfsf8" )
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()








var <- c("Fcer1g", "Lyn", "Syk", "Lat", "Lat2", "Grb2", "Sos1", "Sos2", "Raf-1",
         "Cerk", "Pla2g4a", "Pla2g4b", "Pla2g4c", "Pla2g4f", "Alox5ap", "Ltc4s",
         "Pgd")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis() + labs(title = "Fcer1g pathway -> MAPK pathway" )




