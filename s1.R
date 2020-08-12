### Data Set Processing: h1-h14
### Standardize phenotype data: ha1~ha14
### 2020.08.12
### Song Liting 

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
library(dplyr)
library(stringr)

options(stringsAsFactors = F)
setwd('~/Documents/ScDatabase/processed_data/')

# h1
{
  h1_tpm <- read.table('./1_cortex_development/exprMatrix.txt', sep='\t',header = T,row.names = 1)
  ha1 <- read.table('./1_cortex_development/cotex_m.tsv',sep='\t',header = T,row.names = 1)
  ha1 <- within(ha1,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  
  h1 <- h1_tpm
  
  h1 <- CreateSeuratObject(counts = h1,min.cells  = 5, min.features = 200, meta.data = ha1)
  h1 <- SubsetData(object = h1, cells = rownames(h1@meta.data)[h1@meta.data$Area %in% c("PFC",'MGE','V1')])
  h1[["percent.mt"]] <- PercentageFeatureSet(h1, pattern = "^MT-")
  h1@meta.data$source <- "h1"
  h1@meta.data$tech <- "Fluidigm C1"
  VlnPlot(object = h1, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h1 <- subset(x = h1, subset = nFeature_RNA > 1000 & nFeature_RNA < 6000 )
  
  h1 <- NormalizeData(h1)
  h1 <- FindVariableFeatures(object = h1)
}

# h2
{
  h2_count <- read.table('./2_GSE67835_RAW/h2_count_matrix.txt',sep='\t',header = T)
  rownames(h2_count) <- str_split_fixed(rownames(h2_count),' ',2)[,1]
  ha2 <-   read.table(file='./2_GSE67835_RAW/anno.txt',sep='\t')[,c(1,7,8,9)]
  colnames(ha2) <- c('Sample','Age','Area','Sex')
  ha2$Age_in_Years <- as.numeric(ha2$Age)
  
   ha2 <- within(ha2,{
    Time <- 'P5'
    Time[is.na(Age_in_Years)] <- 'P5'
    Time[Age_in_Years >= 1 & Age_in_Years < 6 ] <- 'P10'
    Time[Age_in_Years >= 6 & Age_in_Years < 12 ] <- 'P11'
    Time[Age_in_Years >= 12 & Age_in_Years < 20 ] <- 'P12'
    Time[Age_in_Years >= 20 & Age_in_Years < 40 ] <- 'P13'
    Time[Age_in_Years >= 40 & Age_in_Years < 60 ] <- 'P14'
    Time[Age_in_Years >= 60  ] <- 'P15'
    
  })
  
  
  h2 <- h2_count
  h2 <- CreateSeuratObject(counts = h2,min.cells  = 5, min.features  = 300, meta.data = ha2,names.field = 1:2, names.delim = "_")
  h2@meta.data$source <- "h2"
  h2@meta.data$tech <- "Fluidigm C1"
  h2[["percent.mt"]] <- PercentageFeatureSet(h2, pattern = "^MT-")
  VlnPlot(object = h2, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)

  h2 <- subset(x = h2, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )

  h2 <- NormalizeData(h2)
  h2 <- FindVariableFeatures(object = h2)
  

}

# h3
{
  h3_count <- read.table('./3_GSE104276_WANG_2017_NATURE/GSE104276_all_pfc_2394_UMI_count_NOERCC.txt',sep='\t',header = T)
  ha3 <- as.data.frame(cbind(colnames(h3_count), substr(colnames(h3_count),1,9), substr(colnames(h3_count),3,4) ))
  colnames(ha3) <- c('cell_name','sample_name','Age_in_Weeks')
  rownames(ha3) <- ha3$cell_name
  ha3$Age_in_Weeks <- as.numeric(ha3$Age_in_Weeks)
  ha3$Area <- 'PFC'
  ha3 <- within(ha3,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  
  haemogolbin <- read.table(file='~/Documents/ScDatabase/processed_data/3_GSE104276_WANG_2017_NATURE/haemogolbin.txt' ,col.names = 'haemogolbin')
  
  h3 <- h3_count
  
  h3 <- CreateSeuratObject(counts = h3,min.cells  = 5, min.features  = 300, meta.data = ha3,names.field = 1:2, names.delim = "_")
  h3 <- SubsetData(h3, cells=names(h3@active.ident)[!names(h3@active.ident) %in% haemogolbin$haemogolbin])
  h3@meta.data$source <- "h3"
  h3@meta.data$tech <- "Smartseq2"
  
  h3[["percent.mt"]] <- PercentageFeatureSet(h3, pattern = "^MT-")
  VlnPlot(object = h3, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h3 <- subset(x = h3, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )

  h3 <- NormalizeData(h3)
  h3 <- FindVariableFeatures(object = h3)
  

}

# h4
{
  h4_tpm <- read.table('./4_GSE103723_WANG_2018_CELLR/h4_tpm_matrix.txt',sep='\t',header = T)
  ha4 <- read.table('./4_GSE103723_WANG_2018_CELLR/anno.txt',sep='\t',header = T)
  abbr_r <- read.table('./4_GSE103723_WANG_2018_CELLR/abbr_region.txt',sep='\t',header = T)
  a4 <- merge(ha4,abbr_r,by='abbr',all.x=T)
  rownames(a4) <- a4$ha4
  ha4 <- within(a4,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  
  h4 <- h4_tpm
  h4 <- CreateSeuratObject(counts = h4,min.cells  = 5, min.features  = 300, meta.data = ha4, names.field = 1:2, names.delim = "_")
  h4@meta.data$source <- "h4"
  h4@meta.data$tech <- "STRT-seq"
  h4[["percent.mt"]] <- PercentageFeatureSet(h4, pattern = "^MT-")
  VlnPlot(object = h4, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h4 <- subset(x = h4, subset = nFeature_RNA > 1000 & nFeature_RNA < 12000 )

  h4 <- NormalizeData(h4)
  h4 <- FindVariableFeatures(object = h4)
  
}

# h5
{
  load('./5_Psychencode_fetal_adult/Sestan.fetalHuman.Psychencode.Rdata')
  h5_count <- count2
  row.names(h5_count) <- str_split_fixed(row.names(h5_count), pattern = '\\|',2)[,2]
  colnames(h5_count) <- str_split_fixed(colnames(h5_count), pattern = '_',2)[,1]
  
  ha5 <- read.table('./5_Psychencode_fetal_adult/scRNA_QC.txt',header = T, sep = '\t',row.names = 2)
  ha5$Age[ ha5$Age=="125days"] <- '18PCW'
  ha5$Age_in_Weeks <- as.numeric(str_split_fixed(ha5$Age,'P',2)[,1])
  
  
  ha5 <- within(ha5,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  
  
  h5 <- h5_count
  h5 <- CreateSeuratObject(counts = h5,min.cells  = 5, min.features  = 300, meta.data = ha5,names.field = 1:2, names.delim = "-")
  h5@meta.data$source <- "h5"
  h5@meta.data$tech <- "Fluidigm C1"
  h5[["percent.mt"]] <- PercentageFeatureSet(h5, pattern = "^MT-")
  
  VlnPlot(object = h5, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h5 <- subset(x = h5, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )

  h5 <- NormalizeData(h5)
  h5 <- FindVariableFeatures(object = h5)
  

}

# h5_adult
{
  load('./5_Psychencode_fetal_adult/Sestan.adultHumanNuclei.Psychencode.Rdata')
  h5_count_a <- umi2
  
  row.names(h5_count_a) <- str_split_fixed(row.names(h5_count_a), pattern = '\\|',2)[,2]

  ha5_a <- read.table('./5_Psychencode_fetal_adult/snRNA_QC.txt',header = T, sep = '\t')
  ha5_a$Cellcode <- tolower(substr(ha5_a$Cellcode,7,8))
  ha5_a$Age_in_Years <- as.numeric(sub('Y','',ha5_a$Age))
  a <- as.data.frame(cbind( colnames(h5_count_a), str_split_fixed(colnames(h5_count_a),'_',2)[,1]))
  
  ha5_a <- merge(ha5_a,a,by.y='V2',by.x='Cellcode')
  
  ha5_a <- within(ha5_a,{
    Time <- 'P5'
    Time[Age_in_Years >= 1 & Age_in_Years < 6 ] <- 'P10'
    Time[Age_in_Years >= 6 & Age_in_Years < 12 ] <- 'P11'
    Time[Age_in_Years >= 12 & Age_in_Years < 20 ] <- 'P12'
    Time[Age_in_Years >= 20 & Age_in_Years < 40 ] <- 'P13'
    Time[Age_in_Years >= 40 & Age_in_Years < 60 ] <- 'P14'
    Time[Age_in_Years >= 60  ] <- 'P15'
  })
  rownames(ha5_a) <- ha5_a$V1
  
  h5_a <- h5_count_a
  h5_a <- CreateSeuratObject(counts = h5_a,min.cells  = 5, min.features  = 300, meta.data = ha5_a)
  h5_a@meta.data$source <- "h5_A"
  h5_a@meta.data$tech <- "10X"
  h5_a[["percent.mt"]] <- PercentageFeatureSet(h5_a, pattern = "^MT-")
  
  VlnPlot(object = h5_a, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h5_a <- subset(x = h5_a, subset = nFeature_RNA > 1000 & nFeature_RNA < 10000 )

  h5_a <- NormalizeData(h5_a)
  h5_a <- FindVariableFeatures(object = h5_a)
  
}

# h7
{
  CerebellarHem <- read.table('./7_GSE97942_RAW/GSE97930_CerebellarHem_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',header = T)
  FrontalCortex <- read.table('./7_GSE97942_RAW/GSE97930_FrontalCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',header = T)
  VisualCortex <- read.table('./7_GSE97942_RAW/GSE97930_VisualCortex_snDrop-seq_UMI_Count_Matrix_08-01-2017.txt',header = T)
  
  h7_count <- transform(Reduce(merge, lapply(list(CerebellarHem,FrontalCortex,VisualCortex), function(x) data.frame(x, rn = row.names(x)))), row.names=rn, rn=NULL)
  colnames(h7_count) <- str_split_fixed(colnames(h7_count),'_',2)[,2]
  h7<- h7_count
  ha7 <- read.table('./7_GSE97942_RAW/anno.txt',sep='\t',header = T,row.names = 5)[,-2]
  ha7 <- within(ha7,{
    Time <- 'P5'
    Time[Age_in_Years >= 1 & Age_in_Years < 6 ] <- 'P10'
    Time[Age_in_Years >= 6 & Age_in_Years < 12 ] <- 'P11'
    Time[Age_in_Years >= 12 & Age_in_Years < 20 ] <- 'P12'
    Time[Age_in_Years >= 20 & Age_in_Years < 40 ] <- 'P13'
    Time[Age_in_Years >= 40 & Age_in_Years < 60 ] <- 'P14'
    Time[Age_in_Years >= 60  ] <- 'P15'
  })
  
  h7 <- CreateSeuratObject(counts = h7,min.cells  = 5, min.features  = 100, meta.data = ha7)
  h7@meta.data$source <- "h7"
  h7@meta.data$tech <- "sndrop-seq"
  h7[["percent.mt"]] <- PercentageFeatureSet(h7, pattern = "^MT-")
  VlnPlot(object = h7, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  dev.print(pdf, file='~/Desktop/h7.pdf')
  
  h7 <- subset(x = h7, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 )

  h7 <- NormalizeData(h7)
  h7 <- FindVariableFeatures(object = h7)
  
}

# h8
{
  load('./8/h8_count.RData')
  ha8 <- read.table('./8/anno.txt',header = T,sep='\t',row.names = 1)
  ha8$Age_in_Weeks <- as.numeric(sub('week_','',ha8$Age_in_Weeks))
  ha8 <- within(ha8,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  
  h8 <- h8_count
  h8 <- CreateSeuratObject(counts = h8,min.cells  = 5, min.features  = 100, meta.data = ha8)
  h8@meta.data$source <- "h8"
  h8@meta.data$tech <- "Fluidigm C1"
  h8[["percent.mt"]] <- PercentageFeatureSet(h8, pattern = "^MT-")
  VlnPlot(object = h8, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h8 <- subset(x = h8, subset = nFeature_RNA > 300 & nFeature_RNA < 5000 )

  h8 <- NormalizeData(h8)
  h8 <- FindVariableFeatures(object = h8)
  
}

# h10
{
  h10_count <- read.table('./10/GTEx_droncseq_hip_pcf/GTEx_droncseq_hip_pcf.umi_counts.txt',header = T,sep='\t')
  h10_count <- h10_count[,grepl('HP|PFC',colnames(h10_count))]
  ha10 <- as.data.frame(colnames(h10_count))
  
  
  ha10 <- as.data.frame(cbind(colnames(h10_count), str_split_fixed(string = colnames(h10_count),pattern = '_|\\.',2)[,1]))
  ha10$Age_in_Years <- '40-65'
  ha10$Time <- 'P14'
  ha10$Sex <- 'Male'
  ha10 <- within(ha10,{
    Area <- 'PFC'
    Area[grepl('HP',V2)] <- 'Hippocampus'
    Area[grepl('PFC',V2)] <- 'PFC'
  })
  rownames(ha10) <- ha10$V1
  
  h10 <- h10_count
  
  h10 <- CreateSeuratObject(counts = h10,min.cells  = 5, min.features  = 100, meta.data = ha10)
  h10@meta.data$source <- "h10"
  h10@meta.data$tech <- "DroNc-seq"
  h10[["percent.mt"]] <- PercentageFeatureSet(h10, pattern = "^MT-")
  
  VlnPlot(object = h10, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0.1)
  h10 <- subset(x = h10, subset = nFeature_RNA > 200 & nFeature_RNA < 3000 )

  h10 <- NormalizeData(h10)
  h10 <- FindVariableFeatures(object = h10)
  
}

# h11
{
  load('/home1/songlt/Documents/ScDatabase/RData/h11_count.RData')
  
  h11 <- h11_d
  h11 <- CreateSeuratObject(counts = h11,min.cells  = 5, min.features  = 100, meta.data = ha11,names.field = 1, names.delim = "_")
  h11@meta.data$source <- "h11"
  h11@meta.data$tech <- "10X"
  h11[["percent.mt"]] <- PercentageFeatureSet(h11, pattern = "^MT-")
  
  VlnPlot(object = h11, features = c("nFeature_RNA", "nCount_RNA"),  ncol = 3, pt.size = 0.1)
  h11 <- subset(x = h11, subset = nFeature_RNA > 300 & nFeature_RNA < 9000 )

  h11 <- NormalizeData(h11)
  h11 <- FindVariableFeatures(object = h11)
  
}

# h12

{
  h12_cpm <- read.table('./12/GSE71315_scell_ncounts.genes.thresh.txt',header = T,row.names = 2)[,-1]
  ha12 <- read.table('./12/anno.txt',header = T)
  rownames(ha12) <- ha12$sample
  ha12 <- within(ha12,{
    Time <- 'P1'
    Time[Age_in_Weeks >= 4 & Age_in_Weeks < 8 ] <- 'P1'
    Time[Age_in_Weeks >= 8 & Age_in_Weeks < 10 ] <- 'P2'
    Time[Age_in_Weeks >= 10 & Age_in_Weeks < 13 ] <- 'P3'
    Time[Age_in_Weeks >= 13 & Age_in_Weeks < 16 ] <- 'P4'
    Time[Age_in_Weeks >= 16 & Age_in_Weeks < 19 ] <- 'P5'
    Time[Age_in_Weeks >= 19 & Age_in_Weeks < 24 ] <- 'P6'
    Time[Age_in_Weeks >= 24 & Age_in_Weeks < 38 ] <- 'P7'
  })
  ha12$Area <- 'NCX'
  
  h12 <- h12_cpm
  h12 <- CreateSeuratObject(counts = h12,min.cells  = 5, min.features  = 100, meta.data = ha12)
  h12@meta.data$source <- "h12"
  h12@meta.data$tech <- "Fluidigm C1"
  h12[["percent.mt"]] <- PercentageFeatureSet(h12, pattern = "^MT-")
  VlnPlot(object = h12, features = c("nFeature_RNA", "nCount_RNA"), ncol = 3, pt.size = 0.1)
  h12 <- subset(x = h12, subset = nFeature_RNA > 300 & nFeature_RNA < 6000 )

  h12 <- NormalizeData(h12)
  h12 <- FindVariableFeatures(object = h12)
  

}

# h13
{
  h13_count <- read.table('./13/GSE81475_Zika.GEO.humanBrain.singleCell.gene.count.txt.gz', sep='\t',header = T, row.names = 1)
  h13_r <- str_split_fixed(rownames(h13_count),'\\|',2)[,2]
  h13_count <- h13_count[!duplicated(h13_r),]
  rownames(h13_count) <- str_split_fixed(rownames(h13_count),'\\|',2)[,2]
  
  ha13 <- read.table('./13/anno.txt',sep='\t')
  
  h13 <- h13_count
  h13 <- CreateSeuratObject(counts = h13,min.cells  = 5, min.features  = 100, meta.data = ha13)
  h13 <- SubsetData(h13, subset.name ='Area', accept.value="Dorsolateral prefrontal cortex")
  
  h13@meta.data$source <- "h13"
  h13@meta.data$tech <- "Fluidigm C1"
  h13[["percent.mt"]] <- PercentageFeatureSet(h13, pattern = "^MT-")
  
  VlnPlot(object = h13, features = c("nFeature_RNA", "nCount_RNA"),  ncol = 3, pt.size = 0.1)
  h13 <- subset(x = h13, subset = nFeature_RNA > 300 & nFeature_RNA < 12000 )

  h13 <- NormalizeData(h13)
  h13 <- FindVariableFeatures(object = h13)
  
}

# h14

{
  ACC <- read.csv('./14/human_ACC_gene_expression_matrices_2018-10-04/human_ACC_2018-10-04_exon-matrix.csv',row.names = 1)
  MTG <- read.csv('./14/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_exon-matrix.csv',row.names = 1)
  VISP <- read.csv('./14/human_VISp_gene_expression_matrices_2018-10-04/human_VISp_2018-10-04_exon-matrix.csv',row.names = 1)
  h14_count <- cbind(ACC, MTG, VISP)
  
  Acc_anno <- read.csv('./14/human_ACC_gene_expression_matrices_2018-10-04/human_ACC_2018-10-04_samples-columns.csv')
  Acc_anno <- Acc_anno[,c('seq_name', "sex", "age_days", "brain_region", "brain_subregion" ,'donor_id')]
  
  MTG_anno <- read.csv('./14/human_MTG_gene_expression_matrices_2018-06-14/human_MTG_2018-06-14_samples-columns.csv')
  MTG_anno <- MTG_anno[,c('sample_name', "sex", "age_days", "brain_region", "brain_subregion" ,'donor_id')]
  colnames(MTG_anno)[1] <- 'seq_name'
  
  VISP_anno <- read.csv('./14/human_VISp_gene_expression_matrices_2018-10-04/human_VISp_2018-10-04_samples-columns.csv')
  VISP_anno <- VISP_anno[,c('seq_name', "sex", "age_days", "brain_region", "brain_subregion" ,'donor_id')]
  
  ha14 <- rbind(Acc_anno,MTG_anno,VISP_anno)
  ha14$Age_in_Years <- as.numeric(ha14$age_days)/365
  rownames(ha14) <- gsub('-','.',ha14$seq_name)
  
  ha14 <- within(ha14,{
    Time <- 'P5'
    Time[Age_in_Years >= 1 & Age_in_Years < 6 ] <- 'P10'
    Time[Age_in_Years >= 6 & Age_in_Years < 12 ] <- 'P11'
    Time[Age_in_Years >= 12 & Age_in_Years < 20 ] <- 'P12'
    Time[Age_in_Years >= 20 & Age_in_Years < 40 ] <- 'P13'
    Time[Age_in_Years >= 40 & Age_in_Years < 60 ] <- 'P14'
    Time[Age_in_Years >= 60  ] <- 'P15'
  })
  
  h14 <- h14_count
  h14 <- CreateSeuratObject(counts = h14,min.cells  = 5, min.features  = 100, meta.data = ha14)
  
  
  h14@meta.data$source <- "h14"
  h14@meta.data$tech <- "smart-seq"
  h14[["percent.mt"]] <- PercentageFeatureSet(h14, pattern = "^MT-")
  
  VlnPlot(object = h14, features = c("nFeature_RNA", "nCount_RNA"),  ncol = 3, pt.size = 0.1)

  h14 <- subset(x = h14, subset = nFeature_RNA > 300 & nFeature_RNA < 10000 )

  h14 <- NormalizeData(h14)
  h14 <- FindVariableFeatures(object = h14)

}

save(h1_tpm,h2_count,h3_count,h4_tpm,h5_tpm,h5_count_a,h7_count,h8_count,h10_count,h11_d,h12_cpm,h13_count,h14_count,file='~/Documents/ScDatabase/RData/h_all_count.RData')
save(ha1,ha2,ha3,ha4,ha5,ha5_a,ha7,ha8,ha10,ha11,ha12,ha13,ha14,file='~/Documents/ScDatabase/RData/h_anno_all.RData')

load(file='~/Documents/ScDatabase/RData/h_anno_all.RData')

# ha1
{
  ha1$priA <- ha1$Area
  ha1$Area[ ha1$Area=='0' & ha1$RegionName=='Cortex'] <- 'Cortex'
  ha1 <- within(ha1,{
    Area[Area=='V1'] <- 'V1C'
  })
  ha1 <- ha1 %>% dplyr::rename( priCluster=WGCNAcluster ,  Sample=Name, Period=Time )
  ha1 <- ha1[,c(-4,-5)]
  ha1$protocal <- 'Fluidigm C1'
  ha1$dataset <- 'h1'
}

#ha2
{
  ha2$priA <- ha2$Area
  ha2$Age_in_Weeks[is.na(ha2$Age_in_Years)] <- 17
  ha2 <- ha2 %>% dplyr::rename( Period=Time )

  ha2 <- within(ha2,{
    Sex[Sex=='M '] <- 'Male'
    Sex[Sex=='F '] <- 'Female'
    Area[Area=='anterior temporal lobe'] <- 'ITC'
    Area[Area=='cortex'] <- 'Cortex'
    
  })
  
  ha2$protocal <- 'Fluidigm C1'
  ha2$dataset <- 'h2'
  
}

#ha3
{
  ha3$priA <- ha3$Area
  ha3 <- ha3[,-1]
  ha3 <- ha3 %>% dplyr::rename( Period=Time ,Sample=sample_name)
  ha3$protocal <- 'Smartseq2'
  ha3$dataset <- 'h3'
}

#ha4 
{
  ha4$Sex <- substr(str_split_fixed( ha4$ha4,'_',4)[,2],4,4)
  ha4 <- within(ha4,{
    Area[Area=='Frontal lobe'] <- 'FC'
    Area[Area=='Temporal lobe'] <- 'TC'
    Area[Area=='Parietal lobe'] <- 'PC'
    Area[Area=='Occipital lobe'] <- 'OC'
    Sex[Sex=='M'] <- 'Male'
    Sex[Sex=='F'] <- 'Female'
  })
  ha4 <- ha4[,c(-1,-2)]
  ha4 <- ha4 %>% dplyr::rename( Period=Time ,Sample=sample, priA=region)
  ha4$protocal <- 'STRT-seq'
  ha4$dataset <- 'h4'
}

# ha5
{
  ha5 <- ha5[,c(2,3,5,8,9,10)]
  
  colnames(ha5) <- c('Sample','Area','Sex','priCluster','Age_in_Weeks','Period')
  ha5_a$priA <- ha5_a$Area
  ha5 <- within(ha5,{
    Area[Area=='pallium'| Area=='Pallium'] <- 'Cortex'
    Sex[Sex=='M'] <- 'Male'
    Sex[Sex=='F'] <- 'Female'
  })
  ha5$protocal <- 'Fluidigm C1'
  ha5$dataset <- 'h5'
}

# ha5_a
{
  ha5_a <- ha5_a[,c(3,4,6,8,10)]
  
  colnames(ha5_a) <- c('Sample','Area','Sex','Age_in_Years','Period')
  ha5_a$priA <- ha5_a$Area
  ha5_a <- within(ha5_a,{
    Sex[Sex=='M'] <- 'Male'
    Sex[Sex=='F'] <- 'Female'
  })
  ha5_a$protocal <- '10X'
  ha5_a$dataset <- 'h5_a'
}

# ha7
{
  ha7 <- ha7[,c(2,3,5:8)]
  
  colnames(ha7) <- c('Age_in_Years','Sex','Sample','Area','priCluster','Period')
  ha7$priA <- ha7$Area
  ha7 <- within(ha7,{
    Area[Area=='Lateral'] <- 'CBC'
    Area[Area=='BA6'] <- 'FC'
    Area[Area=='BA17'] <- 'V1C'
    Area[Area=='BA10'] <- 'FC'
  })
  ha7$protocal <- 'SNdrop-seq'
  ha7$dataset <- 'h7'
}

# ha8
{
  ha8$priA <- ha8$Area
  ha8 <- ha8 %>% dplyr::rename( priCluster=Cell_type ,  Period=Time )
  ha8$Sample <- substr(rownames(ha8),1,7)
  ha8$protocal <- 'Fluidigm C1'
  ha8$dataset <- 'h8'
  ha8$Area <- 'VMB'
}

# ha 10
{
  rownames(ha10) <- ha10$V1
  ha10 <- ha10[,-1]
  colnames(ha10) <- c('Sample','Age_in_Years', 'Period', 'Area')
  ha10$priA <- ha10$Area
  ha10$Area[ha10$Area=='Hippocampus'] <- 'HIP'
  ha10$protocal <- 'DroNc-seq'
  ha10$dataset <- 'h10' 
}

#ha11
{colnames(ha11) <- c('Sample','Sex','Period','Area')
  ha11$priA <- ha11$Area
  ha11$protocal <- '10X'
  ha11$dataset <- 'h11'}

#ha12
{
  ha12$Area <- 'NCX'
  ha12$Sample <- str_split_fixed(ha12$sample,'_|\\.',2)[,1]
  ha12 <- ha12[,-1]
  ha12$protocal <- 'Fluidigm C1'
  ha12$dataset <- 'h12'
  ha12 <- ha12 %>% dplyr::rename(  Period=Time )
  
}

#ha13
{
  ha13$Sample <- substr(rownames(ha13),1,5)
  ha13$priA <- ha13$Area
  ha13$Area[ha13$Area=='Dorsolateral prefrontal cortex'] <- 'DFC'
  ha13 <- ha13 %>% dplyr::rename(  Period=Time )
  ha13$protocal <- 'Fluidigm C1'
  ha13$dataset <- 'h13'
}

#ha14
{
  ha14$Sample <- ha14$donor_id
  
  ha14 <- ha14[,c(2,4,7:9)]
  colnames(ha14) <- c('Sex','Area',"Age_in_Years", "Period",  "Sample"  )
  ha14$priA <- ha14$Area
  ha14 <- within(ha14,{
    Sex[Sex=='M'] <- 'Male'
    Sex[Sex=='F'] <- 'Female'
    Area[Area=='VISp'] <- 'V1C'
  })
  
  ha14$protocal <- 'Smartseq'
  ha14$dataset <- 'h14'
  
}

save(ha1,ha2,ha3,ha4,ha5,ha5_a,ha7,ha8,ha10,ha11,ha12,ha13,ha14,file='~/Documents/ScDatabase/RData/ha_all_norm_0723.RData')
