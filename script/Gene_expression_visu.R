Observation d expresion de gene

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

# 3 Visualisation des genes-----------------------------------------------------
levels(Singlet_norm)


var <- c("Dusp1", "Dusp2" )
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

# Gène unique du venn issu des genes DE aa/VVaa et ab/VVab
var <- c("Dusp1", "Klf4", "Ppp1r15a", "Rgs2" )
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

var <- c("Coq10b" )
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom")

# Gène DE aa & ab et VVaa & VVab
#4 premiers genes aavsab
var <- c("Cd8b1", "Itga4", "Fcer1g", "Ikzf2", "Fcgrt", "Klrk1")
DotPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe") + RotatedAxis()+
  ggtitle("Gènes DE VV aa et VV ab")
VlnPlot(Singlet_norm, features = var, group.by = "pop_Ag_spe", cols = my_palette) + 
  RotatedAxis() +
  theme(legend.position = "bottom") 
