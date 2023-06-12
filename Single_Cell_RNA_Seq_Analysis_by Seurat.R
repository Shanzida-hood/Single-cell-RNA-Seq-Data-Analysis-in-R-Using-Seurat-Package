
## Title :Single cell RNA seq Data Analysis
## Date:05/04/2023
################################################
## about the dataset 
#################################
# script to perform standard workflow steps to analyze single cell RNA-Seq data
# data: 20k Mixture of NSCLC DTCs from 7 donors, 3' v3.1
# data source: https://www.10xgenomics.com/resources/datasets/10-k-human-pbm-cs-multiome-v-1-0-chromium-controller-1-standard-2-0-0         
################################################################################################

setwd("/Users/ssiddi27/Desktop/NGS_data_analysis/Single_Cell_RNA")

# load libraries
library(Seurat)
library(tidyverse)

# Load the NSCLC dataset
nsclc.sparse.m <- Read10X_h5(filename = '/Users/ssiddi27/Desktop/NGS_data_analysis/Single_Cell_RNA/Datasets_single cell/40k_NSCLC_DTC_3p_HT_nextgem_Multiplex_count_raw_feature_bc_matrix.h5')
str(nsclc.sparse.m)
cts <-  nsclc.sparse.m$`Gene Expression`
cts
colnames(cts) ## all barcodes 
rownames(cts) ## genes features/gene name
############################################################################################
# Initialize the Seurat object with the raw (non-normalized data).
nsclc.seurat.obj <- CreateSeuratObject(counts = cts, project = "NSCLC", min.cells = 3, min.features = 200)
str(nsclc.seurat.obj)
nsclc.seurat.obj
# 29552 features across 42081 samples

###########################################################################################

# 1. QC -------
View(nsclc.seurat.obj@meta.data)
# % MT reads
nsclc.seurat.obj[["percent.mt"]] <- PercentageFeatureSet(nsclc.seurat.obj, pattern = "^MT-")
View(nsclc.seurat.obj@meta.data)

VlnPlot(nsclc.seurat.obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(nsclc.seurat.obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") +
  geom_smooth(method = 'lm')
##########################################################################################
# 2. Filtering -----------------
nsclc.seurat.obj <- subset(nsclc.seurat.obj, subset = nFeature_RNA > 200 & nFeature_RNA < 2500 & 
                             percent.mt < 5)
##########################################################################################
# 3. Normalize data ----------
#nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj, normalization.method = "LogNormalize", scale.factor = 10000)
# OR
nsclc.seurat.obj <- NormalizeData(nsclc.seurat.obj)
str(nsclc.seurat.obj)
##########################################################################################

# 4. Identify highly variable features --------------
nsclc.seurat.obj <- FindVariableFeatures(nsclc.seurat.obj, selection.method = "vst", nfeatures = 2000)
##########################################################################################
# Identify the 10 most highly variable genes
top10 <- head(VariableFeatures(nsclc.seurat.obj), 10)
########################################################################################
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(nsclc.seurat.obj)
LabelPoints(plot = plot1, points = top10, repel = TRUE)

#########################################################################################
# 5. Scaling -------------
all.genes <- rownames(nsclc.seurat.obj)
nsclc.seurat.obj <- ScaleData(nsclc.seurat.obj, features = all.genes)

str(nsclc.seurat.obj)
#######################################################################################
# We performed linear dimensionality reductions that is PCA,Principal Component Analysis,to identify the source of heterogenecity in the datasets::
#We run PCA by RunPCA function in provided by seurat object.
## by default it will only consider the variable in seurat object.
# 6. Perform Linear dimensionality reduction --------------
nsclc.seurat.obj <- RunPCA(nsclc.seurat.obj, features = VariableFeatures(object = nsclc.seurat.obj))

# visualize PCA results
print(nsclc.seurat.obj[["pca"]], dims = 1:5, nfeatures = 5)
## we can create a heat map to make sure which principal component shoul we include our downstream analysis by looking at the principal component score and genes that contribute to the principal component.

DimHeatmap(nsclc.seurat.obj, dims = 1, cells = 500, balanced = TRUE)
##in the figure the features show heterogenecity.We will only consider this principal compoent in our downstream analysis.

# determine dimensionality of the data
## that means only choose those statistically significant PCA which capture the majority of the signal in our downstream analysis.This can be done in many way but here we will use elbow plot.
ElbowPlot(nsclc.seurat.obj)
## in the picture x axis represent PC and Y axis represents standard deviation.WE should consider all principal component that explain higher percentage of variance.
## we will consider PCA upto principal component 15.
#######################################################################################
# 7. Clustering ------------

## in this step we will cluster cells that have similar expression pattern.
## for this at first we will find neighbors.
## for this we provide dimensions of principal component first 15.

nsclc.seurat.obj <- FindNeighbors(nsclc.seurat.obj, dims = 1:15)

# understanding resolution
nsclc.seurat.obj <- FindClusters(nsclc.seurat.obj, resolution = c(0.1,0.3, 0.5, 0.7, 1))
View(nsclc.seurat.obj@meta.data)

DimPlot(nsclc.seurat.obj, group.by = "RNA_snn_res.0.5", label = TRUE)

# setting identity of clusters
Idents(nsclc.seurat.obj)
Idents(nsclc.seurat.obj) <- "RNA_snn_res.0.1"
Idents(nsclc.seurat.obj)

# non-linear dimensionality reduction --------------
# If you haven't installed UMAP, you can do so via reticulate::py_install(packages =
# 'umap-learn')
nsclc.seurat.obj <- RunUMAP(nsclc.seurat.obj, dims = 1:15)
# note that you can set `label = TRUE` or use the LabelClusters function to help label
# individual clusters
DimPlot(nsclc.seurat.obj, reduction = "umap")
########################################################################################


