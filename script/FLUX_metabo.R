# Annalyse Flux metabo issu de l'algo porjet long sur data Ly 49 pour test

# 1 Library---------------------------------------------------------------------
if (!require('Seurat')) install.packages('Seurat'); library('Seurat')
library(stringr)
library(openxlsx)


# 2 import data-----------------------------------------------------------------
rm(list = ls())
Singlet_norm <- readRDS("~/Documents/JM/singelcell_LY49/data/Singlet_norm_5.RDS")
predFlux <- read.xlsx('/Users/victor/Documents/JM/singelcell_LY49/data/Ly49_flux_all_module.xlsx', rowNames = T)


all(rownames(predFlux) == colnames(Singlet_norm))

predFlux <- data.matrix(predFlux)
predFlux0 <- t(predFlux)
 
Singlet_norm@assays$FLUX <- CreateAssayObject(counts = predFlux0)  # new groupe dans obj seurat 

Singlet_norm <- RunPCA(Singlet_norm, features = VariableFeatures(object = Singlet_norm), assay = "FLUX")
