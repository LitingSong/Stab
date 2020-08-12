### identification of cell subtypes and marker genes
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
library(RColorBrewer)
library(cowplot)


load('~/Documents/ScDatabase/RData/h.combined1.RData')

# 1. Run the standard workflow for visualization and clustering
DefaultAssay(object = h.combined1) <- "integrated"

h.combined <- ScaleData(object = h.combined1, verbose = FALSE)
h.combined <- RunPCA(object = h.combined, npcs = 40, verbose = FALSE)

## 1.1 t-SNE and Clustering
h.combined <- RunUMAP(object = h.combined, reduction = "pca", dims = 1:30)
h.combined <- RunTSNE(object = h.combined, reduction = "pca", dims = 1:30)
h.combined <- FindNeighbors(object = h.combined, reduction = "pca", dims = 1:30)

h.combined <- FindClusters(h.combined, resolution = 3)
h.combined <- FindClusters(h.combined, resolution = 4)
h.combined <- FindClusters(h.combined, resolution = 2)
DefaultAssay(object = h.combined) <- "RNA"

## 1.2 juege Gender using genes in Y chromosome
meta_data_modify <- function(h.combined){
  
  h.combined@meta.data$USP9Y <- h.combined[['RNA']]@data['USP9Y',]
  h.combined@meta.data$DDX3Y <- h.combined[['RNA']]@data['DDX3Y',]
  h.combined@meta.data$RPS4Y1 <- h.combined[['RNA']]@data['RPS4Y1',]
  meta_data <- h.combined@meta.data
  
  aggregate(meta_data[, c('USP9Y','DDX3Y','RPS4Y1')], list(meta_data$Sample), mean) ->ab
  rownames(ab) <- ab$Group.1
  
  unique(meta_data[,c('Sample','Sex')]) ->sex
  rownames(sex) <- sex$Sample
  
  ab$sex <- sex[rownames(ab),'Sex']
  # 设置阈值MAle 这三个基因中至少有一个表达>0.2,Female这三个基因表达全部<0.2
  ab <- within(ab, {
    Sex <- 'Male'
    Sex[USP9Y <0.2 & DDX3Y<0.2 & RPS4Y1<0.2] <- 'Female'
    Sex[USP9Y >0.2 | DDX3Y>0.2 | RPS4Y1>0.2] <- 'Male'
  })
  
  # 设置阈值MAle 这三个基因中至少有一个表达>0.2,Female这三个基因表达全部<0.2
  h.combined@meta.data <- within(h.combined@meta.data, {
    Sex <- 'Male'
    Sex[Sample%in% rownames(ab)[ab$Sex=='Female'] ] <- 'Female'
    Sex[Sample%in% rownames(ab)[ab$Sex=='Male'] ] <- 'Male'
  })
  
  
  # 年龄全部转为天
  h.combined@meta.data <- within(h.combined@meta.data, {
    Age_in_days <- 0
    Age_in_days[!is.na(Age_in_Weeks)] <- 7*as.numeric(Age_in_Weeks[!is.na(Age_in_Weeks)])
    Age_in_days[!is.na(Age_in_Years)] <- 365*as.numeric(Age_in_Years[!is.na(Age_in_Years)])
  })
  
  return(h.combined)
}
h.combined <- meta_data_modify(h.combined=h.combined)

cls <- sort(unique(h.combined$integrated_snn_res.2))

# 2. define cell type based on marker genes

df <- data.frame(row.names = cls)
for(gene in c('SLC17A7','GAD1','HES1','PDGFRA','PCDH15','GFAP','SLC1A2','MBP','PTPRC','FLT1','PDGFRB')){
  for(cl in cls){
    df[as.character(cl), gene] <- mean(subset( h.combined, integrated_snn_res.2==cl)@assays$RNA[gene,])
  }
}

for(cl in cls){
  df[as.character(cl), 'cl' ] <- c('ExN','InN','NPC','OPC','OPC','Astro','Astro','Olig','Micro','Endo','Perc' )[which.max(df[as.character(cl),])]
}

load('~/Documents/ScDatabase/RData/df_res2_pca30.RData')

o.ids <- as.character(cls)
n.ids <- paste(df$cl, cls)
h.combined@active.ident <- plyr::mapvalues(x = h.combined@active.ident, from = o.ids, to = n.ids)

# 3. remove cell subtypes with higher expression in both neuronal and non-neuronal cells
'%ni%' <- Negate('%in%')
h.combined <- subset(h.combined, seurat_clusters%ni%c('53','20'))

## 4.1 merge cell stubtypes 1
o.ids <- c(c('Olig 0',  'Olig 1', 'Olig 25', 'Olig 51'), c('Olig 5', 'Olig 24', 'Olig 19'), 'Astro 4','Astro 9', 'ExN 21','InN 39')
n.ids <- c(rep('Olig_0_1_25_51',4),rep('Olig_5_24_19',3), 'Astro_4_9','Astro_4_9','Gran','Purk')
h.combined@active.ident <- plyr::mapvalues(x = h.combined@active.ident, from = o.ids, to = n.ids)

## find marker
DefaultAssay(object = h.combined) <- "RNA"
h.combined_ExN <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('ExN',h.combined@active.ident)])
h.combined_InN <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('InN',h.combined@active.ident)])
h.combined_Astro <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('Astro',h.combined@active.ident)])
h.combined_Olig <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('Olig',h.combined@active.ident)])
h.combined_opc <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('OPC',h.combined@active.ident)])

h.combined_ExN.markers <- FindAllMarkers(object = h.combined_ExN,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_InN.markers <- FindAllMarkers(object = h.combined_InN,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_Astro.markers <- FindAllMarkers(object = h.combined_Astro,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_Olig.markers <- FindAllMarkers(object = h.combined_Olig,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_opc.markers <- FindAllMarkers(object = h.combined_opc,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)


## 4.2. merge2
for(mf in list(h.combined_Astro.markers, h.combined_ExN.markers,h.combined_InN.markers,h.combined_Olig.markers,h.combined_opc.markers)){
  mf <- subset(mf,avg_logFC>0.69 & p_val_adj<0.01 )
  mf$cluster <- as.character(mf$cluster)
  marker_freq <- as.data.frame(table(mf$gene))
  marker_freq <- merge(mf, marker_freq, by.x='gene',by.y='Var1',all.x=T) 
  speci_marker <- subset(marker_freq, Freq==1)
  print(table(speci_marker$cluster))
}

h.combined_tree <- BuildClusterTree(object = h.combined, dims = 1:30)
PlotClusterTree(object = h.combined_tree)
DimPlot(object = h.combined, reduction = "umap", label=T)+NoLegend()

o.ids <- c(c('Olig_0_1_25_51',  'Olig 47', 'Olig 49'), c( 'ExN 18','ExN 50', 'ExN 2') ,'ExN 43','ExN 14')
n.ids <- c(rep('Olig_0_1_25_51_47_49',3),rep('ExN_18_50_2',3), rep('ExN_43_14',2))
h.combined@active.ident <- plyr::mapvalues(x = h.combined@active.ident, from = o.ids, to = n.ids)

h.combined_ExN <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('ExN',h.combined@active.ident)])
h.combined_Olig <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('Olig',h.combined@active.ident)])

h.combined_ExN.markers <- FindAllMarkers(object = h.combined_ExN,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_Olig.markers <- FindAllMarkers(object = h.combined_Olig,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)

for(mf in list( h.combined_ExN.markers,h.combined_Olig.markers)){
  mf <- subset(mf,avg_logFC>0.69 & p_val_adj<0.01 )
  mf$cluster <- as.character(mf$cluster)
  marker_freq <- as.data.frame(table(mf$gene))
  marker_freq <- merge(mf, marker_freq, by.x='gene',by.y='Var1',all.x=T) 
  speci_marker <- subset(marker_freq, Freq==1)
  print(table(speci_marker$cluster))
}

# 5. identification of marker genes
DefaultAssay(object = h.combined) <- "RNA"

h.combined_ExN <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('ExN',h.combined@active.ident)])
h.combined_InN <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('InN',h.combined@active.ident)])
h.combined_Astro <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('Astro',h.combined@active.ident)])
h.combined_Olig <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('Olig',h.combined@active.ident)])
h.combined_opc <- SubsetData(h.combined, cells=colnames(h.combined)[grepl('OPC',h.combined@active.ident)])

h.combined_ExN.markers <- FindAllMarkers(object = h.combined_ExN,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_InN.markers <- FindAllMarkers(object = h.combined_InN,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_Astro.markers <- FindAllMarkers(object = h.combined_Astro,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_Olig.markers <- FindAllMarkers(object = h.combined_Olig,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
h.combined_opc.markers <- FindAllMarkers(object = h.combined_opc,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)

h.combined.markers <- FindAllMarkers(object = h.combined,  thresh.use = 0.25, only.pos=T, logfc.threshold = 0.4)
save(h.combined.markers, h.combined_ExN.markers, h.combined_InN.markers ,h.combined_Astro.markers ,h.combined_Olig.markers,h.combined_opc.markers ,file='/home1/songlt/Documents/ScDatabase/RData/markers_res2_pc30_refine.RData' )


save(h.combined, file='~/Documents/ScDatabase/RData/h.combined_pc30_refine.RData')
save(h.combined.markers, h.combined_ExN.markers, h.combined_InN.markers ,h.combined_Astro.markers ,h.combined_Olig.markers,h.combined_opc.markers ,file='/home1/songlt/Documents/ScDatabase/RData/markers_res2_pc30_refine.RData' )
load('~/Documents/ScDatabase/RData/h.combined_pc30_refine.RData')

h.combined <- subset(h.combined, seurat_clusters%ni%c('53','20','48'))

o.ids <- sort(levels(h.combined@active.ident))

n.ids <- c(paste('Astro',1:4,sep=''), 'Endo', paste('ExN',1:15,sep=''), "Gran", paste('InN',1:10,sep=''), "Micro",'NPC',paste('Olig',1:4,sep=''),paste('OPC',1:3,sep=''), "Perc", "Purk" )
as.data.frame(cbind(o.ids,n.ids))->ids
rownames(ids) <- ids$o.ids
save(ids,file='~/Desktop/ids.RData')

h.combined@active.ident <- plyr::mapvalues(x = h.combined@active.ident, from = o.ids, to = n.ids)
h.combined$cluster1 <- h.combined@active.ident
save(h.combined, file='~/Documents/ScDatabase/RData/h.combined_pc30_refine2.RData')

load(file='~/Documents/ScDatabase/RData2/h.combined_pc30_refine2.RData') 


# 6. annotation of cell subtypes
lake_marker <- read.table('~/Documents/ScDatabase/processed_data/lake_2015_marker.txt',sep='\t',header = T)

## 6.1. INn 
h.combined$cluster <- gsub('\\d','',h.combined$cluster1)

lake_in <-  lake_marker[lake_marker$cluster %in%c(paste('In', 1:8,sep='')) ,'Gene']

#h.combined_InN <- subset(h.combined,cluster=='InN' & dataset =='h7')
h.combined_InN <- subset(h.combined,cluster=='InN' & Period%in%c('P12','P13','P14','P15'))
h.combined_InN <- subset(h.combined,cluster=='InN' )


myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

myplot <- DotPlot(h.combined_InN, features = rev(lake_in),cols=  c('khaki1',myPalette(100))) + RotatedAxis() +
  theme(legend.position='top' , axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.75)))

ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/FigS5_inn_lake.png', 
       plot=myplot,width = 8.5, height = 4.5,units = 'in',dpi=300)

## 6.2 Exn 
lake_ex <-  unique(lake_marker[lake_marker$cluster %in%c(paste('Ex', 1:8,sep='')) ,'Gene'])
h.combined_ExN <- subset(h.combined,cluster=='ExN' &  Period%in%c('P12','P13','P14','P15'))

h.combined_ExN@active.ident <- factor(h.combined_ExN@active.ident, levels = paste('ExN',1:15,sep=''))

myplot <- DotPlot(h.combined_ExN, features = rev(lake_ex),cols=  c('khaki1',myPalette(100))) + RotatedAxis() +  
  theme(legend.position='' , axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.3)))


ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/FigS5_exn_lake.png', 
       plot=myplot,width = 8.5, height = 5.5,units = 'in',dpi=300)


o.ids <- c(paste('InN', 1:10, sep=''), paste('ExN', 1:15,sep=''))
n.ids <- c(paste('InN', c(3, 6, 5, '1a', '4a', '4_5', '1b', '4b', '5_6', '7_8' ), sep=''),
           paste('ExN', c('1_4', '1a', 3, '6a', '2', '9', 10, '6b', 11, '8', '4_6', '1b', 5, '1c', 4 ), sep=''))

h.combined@active.ident <- plyr::mapvalues(x = h.combined@active.ident, from = o.ids, to = n.ids)


Area_levels <- c('NCX' ,'FC','DFC','PFC', 'OC', 'V1C' , 'PC' ,'TC', 'ITC', 'MTG', 'CTX', 'MGE', 'ACC', 'VMB','HIP','CBC','IG','SN', 'MDL','Pons')
cluster_levels <- c('NPC', paste('ExN',c("1a" , "1b" , "1c", "1_4", "2"   ,"3"  , "4" ,  "4_6" ,"5" ,  "6a",  "6b" , "8" ,  "9", "10",  "11"  ),sep=''),
                    paste('InN', c("1a",  "1b" , "3" ,  "4a" , "4b",  "4_5", "5",   "5_6" ,"6",  "7_8"), sep=''),
                    paste('Astro',1:length(grep('Astro',levels(h.combined$cluster1))),sep=''),
                    paste('Olig',1:length(grep('Olig',levels(h.combined$cluster1))),sep=''),
                    paste('OPC',1:length(grep('OPC',levels(h.combined$cluster1))),sep=''),
                    'Micro','Endo','Perc','Gran','Purk')

h.combined@active.ident <- factor(h.combined@active.ident, levels = cluster_levels)
h.combined$cluster <- gsub('\\d','',h.combined$cluster1)
h.combined$cluster1 <- h.combined@active.ident
setwd('~/Documents/ScDatabase/processed_data/')

cell_type_info <- read.table('./cell_types.txt',row.names = 1, sep='\t',header = T)
period_infos <- read.table('./Period.txt',header = T,sep='\t',row.names = 1)
Area_infos <- read.table('./Area.txt',header = T,sep='\t',row.names = 1)

h.combined$cluster1 <- factor(h.combined$cluster1,levels=cluster_levels)
h.combined$Period <- factor(h.combined$Period,levels=c(paste('P',1:15,sep='')))
h.combined$dataset <- factor(h.combined$dataset,levels=c(paste('h',1:5,sep=''),'h5_a',paste('h',c(7,8,10:14),sep='') ))
h.combined$Area <- factor(h.combined$Area,levels =Area_levels )
cls <- c("NPC", "ExN" , "InN" , "OPC" , "Olig", "Astro",  "Micro","Endo",  "Perc",   "Gran" ,       "Purk"  )
h.combined$cluster <- factor(h.combined$cluster,levels =cls ) 

save(h.combined, file='~/Documents/ScDatabase/RData/h.combined_pc30_refine3.RData')


load('~/Documents/ScDatabase/RData/markers_res2_pc30_refine22.RData')
h.combined.markers$cluster <- plyr::mapvalues(x = h.combined.markers$cluster, from = o.ids, to = n.ids)
h.combined.markers$cluster <- factor(h.combined.markers$cluster, levels=cluster_levels)
save(h.combined.markers, file='~/Documents/ScDatabase/RData/markers_res2_pc30_refine3.RData')
load('~/Documents/ScDatabase/RData/markers_res2_pc30_refine3.RData')

h.combined.markers <- subset(h.combined.markers, p_val_adj<0.05 )
marker_freq <- as.data.frame(table(h.combined.markers$gene))
h.combined.markers <- merge(h.combined.markers, marker_freq, by.x='gene',by.y='Var1',all.x=T) 
h.combined.markers <- subset(h.combined.markers, Freq <=4)
save(h.combined.markers, file='~/Documents/ScDatabase/RData/markers_res2_pc30_refine4.RData')

cell_type_info <- read.table('./cell_types.txt',row.names = 1, sep='\t',header = T,stringsAsFactors = F)
o.ids <- c(paste('InN', 1:10, sep=''), paste('ExN', 1:15,sep=''))
n.ids <- c(paste('InN', c(3, 6, 5, '1a', '4a', '4_5', '1b', '4b', '5_6', '7_8' ), sep=''),
           paste('ExN', c('1_4', '1a', 3, '6a', '2', '9', 10, '6b', 11, '8', '4_6', '1b', 5, '1c', 4 ), sep=''))
rownames(cell_type_info) <- plyr::mapvalues(x = rownames(cell_type_info), from = o.ids, to = n.ids)
cell_type_info[n.ids,'Name'] <- n.ids
cell_type_info[, 'Name'] <- gsub('ExN','Excitatory Neuron' ,cell_type_info[, 'Name'])
cell_type_info[, 'Name'] <- gsub('InN','InterNeuron' ,cell_type_info[, 'Name'])
write.table(cell_type_info, file='~/Documents/ScDatabase/processed_data/cell_types_v2.txt', sep='\t')


