### Generate seurat objects for each dataset
### combiend data sets
### remove batch effects
### 2020.08.12
### Song Liting 

library(scran)
library(Seurat)
library(mclust)
library(dplyr)
library(SingleCellExperiment)
library(sctransform)
library(scater)
library(stringr)
library(plotly)
library(dbscan)
library(reshape2)
library(Matrix)
library(pheatmap)

options(stringsAsFactors = F)
setwd('~/Documents/ScDatabase/processed_data/')
load('~/Documents/ScDatabase/RData/ha_all_norm_0723.RData')
load('~/Documents/ScDatabase/RData/h_all_count.RData')


gene_length <- read.table('./gene_length.txt',header = T,row.names = 2)
genes <- rownames(gene_length)

# 1. Generate seurat objects for each dataset
## h1
{
  h1 <- CreateSeuratObject(counts = h1_tpm[intersect(genes,rownames(h1_tpm)),],min.cells  = 5, min.features = 200, meta.data = ha1)
  h1 <- SubsetData(object = h1, cells = rownames(h1@meta.data)[h1@meta.data$Area %in% c("PFC",'MGE','V1')])
  h1 <- subset(x = h1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 )
  h1 <- NormalizeData(h1)
  h1 <- FindVariableFeatures(object = h1)
}

## h2
{
  
  h2 <- CreateSeuratObject(counts = h2_count[intersect(genes,rownames(h2_count)),], min.cells  = 5, min.features  = 300, meta.data = ha2,names.field = 1:2, names.delim = "_")
  h2 <- subset(x = h2, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )
  h2 <- NormalizeData(h2)
  h2 <- FindVariableFeatures(object = h2)

}

## h3
{
  haemogolbin <- read.table(file='~/Documents/ScDatabase/processed_data/3_GSE104276_WANG_2017_NATURE/haemogolbin.txt' ,col.names = 'haemogolbin')
  h3 <- CreateSeuratObject(counts = h3_count[intersect(genes,rownames(h3_count)),], min.cells  = 5, min.features  = 300, meta.data = ha3,names.field = 1:2, names.delim = "_")
  h3 <- SubsetData(h3, cells=names(h3@active.ident)[!names(h3@active.ident) %in% haemogolbin$haemogolbin])
  h3 <- subset(x = h3, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )
  h3 <- NormalizeData(h3)
  h3 <- FindVariableFeatures(object = h3)

}

## h4
{
  h4 <- CreateSeuratObject(counts = h4_tpm[intersect(genes,rownames(h4_tpm)),],min.cells  = 5, min.features  = 300, meta.data = ha4, names.field = 1:2, names.delim = "_")
  h4 <- subset(x = h4, subset = nFeature_RNA > 1000 & nFeature_RNA < 12000 )
  h4 <- NormalizeData(h4)
  h4 <- FindVariableFeatures(object = h4)
}

## h5
{
  h5 <- CreateSeuratObject(counts = h5_tpm[intersect(genes,rownames(h5_tpm)),],min.cells  = 5, min.features  = 300, meta.data = ha5,names.field = 1:2, names.delim = "-")
  h5 <- subset(x = h5, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )
  h5 <- NormalizeData(h5)
  h5 <- FindVariableFeatures(object = h5)

}

## h5_adult
{
  
  h5_a <- CreateSeuratObject(counts = h5_count_a[intersect(genes,rownames(h5_count_a)),],min.cells  = 5, min.features  = 300, meta.data = ha5_a)
  h5_a <- subset(x = h5_a, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )
  h5_a <- NormalizeData(h5_a)
  h5_a <- FindVariableFeatures(object = h5_a)
}

## h7
{
  
  h7 <- CreateSeuratObject(counts = h7_count[intersect(genes,rownames(h7_count)),],min.cells  = 5, min.features  = 100, meta.data = ha7)
  h7 <- subset(x = h7, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 )
  h7 <- NormalizeData(h7)
  h7 <- FindVariableFeatures(object = h7)
}

## h8
{
  h8 <- CreateSeuratObject(counts = h8_count[intersect(genes,rownames(h8_count)),],min.cells  = 5, min.features  = 100, meta.data = ha8)
  h8 <- subset(x = h8, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 )
  h8 <- NormalizeData(h8)
  h8 <- FindVariableFeatures(object = h8)
}

## h10
{
  h10 <- CreateSeuratObject(counts = h10_count[intersect(genes,rownames(h10_count)),],min.cells  = 5, min.features  = 100, meta.data = ha10)
  h10 <- subset(x = h10, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 )
  h10 <- NormalizeData(h10)
  h10 <- FindVariableFeatures(object = h10)
}

## h11
{
  
  h11 <- CreateSeuratObject(counts = h11_d[intersect(genes,rownames(h11_d)),],min.cells  = 5, min.features  = 100, meta.data = ha11,names.field = 1, names.delim = "_")
  h11 <- subset(x = h11, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 )
  h11 <- NormalizeData(h11)
  h11 <- FindVariableFeatures(object = h11)

}

## h12
{
  h12 <- CreateSeuratObject(counts = h12_cpm[intersect(genes,rownames(h12_cpm)),],min.cells  = 5, min.features  = 100, meta.data = ha12)
  h12 <- subset(x = h12, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 )
  h12 <- NormalizeData(h12)
  h12 <- FindVariableFeatures(object = h12)
}

## h13
{
  
  h13 <- CreateSeuratObject(counts = h13_count[intersect(genes,rownames(h13_count)),], min.cells  = 5, min.features  = 100, meta.data = ha13)
  h13 <- SubsetData(h13, subset.name ='Area', accept.value="DFC")
  h13 <- subset(x = h13, subset = nFeature_RNA > 300 & nFeature_RNA < 12000 )
  h13 <- NormalizeData(h13)
  h13 <- FindVariableFeatures(object = h13)
}

## h14
{
  h14 <- CreateSeuratObject(counts = h14_count[intersect(genes,rownames(h14_count)),], min.cells  = 5, min.features  = 100, meta.data = ha14)
  h14 <- subset(x = h14, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 )
  h14 <- NormalizeData(h14)
  h14 <- FindVariableFeatures(object = h14)

}

save(h1,h2,h3,h4,h5,h5_a,h7,h8,h10,h11,h12,h13,h14,file='~/Documents/ScDatabase/RData/h_seurat_v0723.RData')

# 2. integrated data and remove batch effect

h.anchors <- FindIntegrationAnchors(object.list = list(h1,h2,h3,h4,h5,h5_a,h7,h8,h10,h11,h12,h13,h14))
save(h.anchors,file='~/Documents/ScDatabase/RData/h.anchors_v0723.RData')
h.combined1 <- IntegrateData(anchorset = h.anchors, dims = 1:40)
save(h.combined1,file='~/Documents/ScDatabase/RData/h.combined1.RData')

