## Figuers in munuscript

library(Seurat)
library(reshape2)
library(dplyr)
library(RColorBrewer)
library(cowplot)
library(base2grob)
library(gridExtra)
library(stringr)

options(stringsAsFactors = F)

load('~/Documents/ScDatabase/RData/h.combined_pc30_refine3.RData')
meta_data <- h.combined@meta.data

setwd('~/Documents/ScDatabase/processed_data/')
cell_type_info <- read.table('./cell_types_v2.txt',row.names = 1, sep='\t',header = T)
period_infos <- read.table('./Period.txt',header = T,sep='\t',row.names = 1)
Area_infos <- read.table('./Area.txt',header = T,sep='\t',row.names = 1)
o.ids <- c(paste('InN', 1:10, sep=''), paste('ExN', 1:15,sep=''))
n.ids <- c(paste('InN', c(3, 6, 5, '1a', '4a', '4_5', '1b', '4b', '5_6', '7_8' ), sep=''),
           paste('ExN', c('1_4', '1a', 3, '6a', '2', '9', 10, '6b', 11, '8', '4_6', '1b', 5, '1c', 4 ), sep=''))
meta_d <- as.data.frame(table(meta_data[,c('Area','Period','cluster1')]))
meta_d <- meta_d[meta_d$Freq!=0,]


#############
## Table 1
#############

Table1 <- table(meta_data$dataset)

ha <- unique(meta_data[, c('Area','Period','dataset','protocal')])

aa <- aggregate(ha$Area, list(ha$dataset), paste, collapse=", ") 

dataset_sum <- cbind(as.data.frame(table(meta_data$dataset)),
                     as.data.frame(ha %>%group_by(dataset) %>%summarise(Area=paste(unique(Area),collapse=','))),
                     as.data.frame(ha %>%group_by(dataset) %>%summarise(Period=paste(unique(Period),collapse=','))),
                     as.data.frame(ha %>%group_by(dataset) %>%summarise(Protocol=paste(unique(protocal),collapse=','))))[,c(1,2,4,6,8)]

dataset_sum1 <- read.table('~/Documents/ScDatabase/Tables/sample.txt',sep='\t',comment.char = '',header = T)
dataset_sum1[,c('X..cells','Areas')] <- dataset_sum[,c('Freq','Area')]
write.table(dataset_sum1, file='~/Documents/ScDatabase/Tables/Table1_data_set_sum.txt', sep='\t',quote = F, row.names = F)

#############
## Figure 1
#############

# Umap
p1 <- DimPlot(h.combined)+NoLegend() +theme_bw()+theme(axis.text.x = element_blank(),
                                                       axis.text.y=element_blank(),
                                                       axis.ticks.y =element_blank(),
                                                       axis.ticks.x =element_blank(),
                                                       axis.title.x=element_blank(),
                                                       axis.title.y=element_blank(),
                                                       legend.position="none",
                                                       panel.background=element_blank(),
                                                       panel.border=element_blank(),
                                                       panel.grid.major=element_blank(),
                                                       panel.grid.minor=element_blank(),
                                                       plot.background=element_blank())

# marker
plot_vln <- function(t, gene) {
  d <- as.matrix(t[['RNA']]@data[intersect(gene, rownames(t[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$tree.ident <- t$cluster1[dd$cell]
  #str(dd$tree.ident)
  dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  
  ggplot(dd, aes(tree.ident, value, fill = tree.ident)) + geom_violin(scale = "width", trim = T, alpha = 0.8, adjust = 1) + 
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
    theme(strip.background = element_blank(), strip.placement = "outside", 
          strip.text.y = element_text(size = 10,face ='italic'),  panel.grid = element_blank(), axis.title.y = element_blank(),
          panel.border = element_blank()) + 
    theme(strip.text.y = element_text( size = 10),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = rel(0.9)),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    scale_y_continuous(limits=c(0,5))
}

p2 <- plot_vln(t=subset(h.combined, cluster1%in%c('ExN1a','ExN1b','InN1a','InN1b')), gene=c( 'SLC17A7', 'GAD1') )

# pie
p3 <-  base2grob(~pie(c(2,2,4,1),labels ='',col=c('yellow3','tan1','yellowgreen','peachpuff')))
# composition P3
meta_d <- as.data.frame(table(meta_data[,c('Area','Period','cluster1')]))
meta_d <- meta_d[meta_d$Freq!=0,]
CEll_Types_Periods <- function(period){

  as<- meta_d[meta_d$Period==period,]
  as <- subset(as, cluster1%in%c('NPC','ExN1b','InN1b'))
  
  as <- as %>% group_by(Area) %>%
    dplyr::mutate(percent = Freq/sum(Freq)) %>% as.data.frame()
  

  
  p4 <- ggplot(as, aes(x = (cluster1), y = Area)) +
    geom_point(aes(size= percent) ,color= 'yellowgreen' ) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90,hjust=1),legend.position = '')+
    scale_color_manual(values = c('yellowgreen','tan1'))+
    ylab("") +
    xlab("") +
    ggtitle(paste(period,' (',period_infos[period,'age'], ')' ,sep=''))
  
}
CEll_Types_Periods(period='P3')

# expression P3 
Period_geneE <- function(period, gene){
  
  h.combinedP3 <- subset( h.combined, Period==period & cluster1%in%c('NPC','InN5') )
  d <- as.matrix(h.combinedP3[['RNA']]@data[gene,])
  dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
  dd$Gene <- gene
  dd$cell <- as.character(dd$cell)
  dd$cluster <- meta_data[dd$cell,'cluster1']
  
  dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
  dd<- dd[!is.na(dd$cluster),]
  dd <- dd[dd$Period==period,]
  
  p5 <- ggplot(dd,aes(x=Area, y=value, fill=Area, color=Area,alpha=0.4)) +
    geom_violin( size=0.5) + 
    geom_boxplot(width=0.1, outlier.size = 0.1)+
    facet_wrap(~cluster) +
    theme_light()+
    theme(legend.position='' , axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = rel(0.9)))+
    ylab('Expression Value') +
    xlab('') +
    ggtitle(paste(gene, ' (', period,')',collapse = '',sep=''))
 
}
Period_geneE(gene='GAD1',period='P3')

# composition FC
CEll_Types_Area <- function(area){
  
  as<- meta_d[meta_d$Area==area,]
  as <- subset(as, cluster1%in%c('NPC','ExN1a','InN1a'))
  
  as <- as %>% group_by(Period) %>%
    dplyr::mutate(Proportion = Freq/sum(Freq)) %>% as.data.frame()
  
  p6 <- ggplot(as, aes(x = (cluster1), y = Period)) +
    geom_point(aes(size= Proportion ),color='yellowgreen') +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5),legend.position = '')+
    scale_color_manual(values = c('yellowgreen','tan1'))+
    ylab("") +
    xlab("") +
    ggtitle(paste(area,' (',Area_infos[area,'Area'], ')' ,sep=''))
  
}
CEll_Types_Area(area='FC')

# expression FC
Area_geneE <- function(Region, gene){
  
  #d <- as.matrix(h.combined[['RNA']]@data[gene,])
  h.combinedFC <- subset( h.combined, Area==Region & cluster1%in%c('ExN1_4','InN5','InN7_8') )
  
  d <- as.matrix(h.combinedFC[['RNA']]@data[gene,])
  
  dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
  dd$Gene <- gene
  dd$cell <- as.character(dd$cell)
  dd$cluster <- meta_data[dd$cell,'cluster1']
  dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
  dd<- dd[!is.na(dd$cluster),]
  dd <- dd[dd$Area==Region,]
  
  p7 <- ggplot(dd,aes(x=Period, y=value, color=cluster,group=cluster,fill=cluster)) +
    geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15)+
    xlab('Period') +
    ylab('Expression Value')+
    ggtitle(paste(Region," (", gene, ")", collapse = '',sep=''))+
    theme_light()+
    theme(legend.key.size = unit(0.5, "cm")) 
  


}
Area_geneE(gene='GAD1',Region='FC')

#############
### Figure 1
#############


P6_Seurat <- subset( h.combined, Period=='P6' & Area=='PC')
p6_df <- as.data.frame(P6_Seurat[['RNA']]@data[c('SOX2',    'SLC17A7', 'GAD1',  'AQP4' , 'MBP'  ,   'PCDH15', 'APBB1IP', 'FLT1'  , 'PDGFRB'  ) ,]) 
anno_col_all <- P6_Seurat@meta.data[,c('cluster1','cluster')]

p8 <- pheatmap(p6_df[,rownames(anno_col_all)],cluster_rows=F,cluster_cols=F,show_colnames=F,show_rownames=T,fontsize_row=7,annotation_col=anno_col_all)


#############
## Figure 2
#############
myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

myplot <- DimPlot(h.combined,label = T,label.size = 2)+NoLegend()

cell_type_info <- read.table('./cell_types.txt',row.names = 1, sep='\t',header = T)

plot_vln <- function(t, gene) {
  d <- as.matrix(t[['RNA']]@data[intersect(gene, rownames(t[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$tree.ident <- t$cluster1[dd$cell]
  #str(dd$tree.ident)
  dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  
  ggplot(dd, aes(tree.ident, value, fill = tree.ident)) + geom_violin(scale = "width", trim = T, alpha = 0.8, adjust = 1) + 
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "left") +
    theme(strip.background = element_blank(), strip.placement = "outside", 
          strip.text.y = element_text(colour = "black", size = 8, face ='italic'),  panel.grid = element_blank(), axis.title.y = element_blank(),
          panel.border = element_blank()) + 
    theme(strip.text.y = element_text(size = 8,angle = 180),
          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    scale_y_continuous(limits=c(0,5))
}

myplot <- plot_vln(t=h.combined, gene=c('SOX2',    'SLC17A7', 'GAD1',  'AQP4' , 'MBP'  ,   'PCDH15', 'APBB1IP', 'FLT1'  , 'PDGFRB'  ) )


layer_marker <- read.table('~/Documents/ScDatabase/Tables/layer_marker1.txt',sep='\t',header = T,row.names = 1)
layer_marker <- layer_marker[layer_marker$Cortical.marker..human. !='',]
layer_marker <- layer_marker[ rownames(layer_marker)%in%rownames(h.combined@assays$RNA) , ]

layer_marker <- layer_marker[ rownames(layer_marker)%in%rownames(h.combined@assays$RNA) , ]
layer_marker <- layer_marker[order(layer_marker$`Cortical.marker..human.` ), ]

layer_marker <- layer_marker[grepl('layer', layer_marker$`Cortical.marker..human.` ), ]
layer_marker$L <- gsub( 'interneuron/','',layer_marker$`Cortical.marker..human.` )
layer_marker$L <- gsub( 'layer|layer\\s','L', layer_marker$L)
layer_marker$L <- gsub( '/astrocyte','', layer_marker$L)

## ExN
h.combined_ExN <- subset(h.combined,cluster=='ExN')
ExN_layer <- c('CNR1','LAMP5','CBLN2','CUX2','RORB','TOX','ETV1','TLE4','NR4A2','NTNG2')

myplot <- DotPlot(h.combined_ExN, features = (ExN_layer) )+ RotatedAxis() +coord_flip() +
  scale_x_discrete(labels = rev(paste( layer_marker[(ExN_layer),'L'], (ExN_layer))) )+xlab('')+ylab('')+NoLegend()+
  theme(axis.text.y = element_text(size = 8,face ='italic'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

## InN
InN_layer <- (c("CXCL14", "CHRNA7" ,'CNR1',   "LAMP5", "SV2C" ,  'SNCG', "SULF1" ,'TOX','PDE1A', 'SYNPR'))
h.combined_InN <- subset(h.combined,cluster=='InN')
myplot <- DotPlot(h.combined_InN, features = (InN_layer)) + RotatedAxis() +coord_flip() + 
  scale_x_discrete(labels = rev(paste( layer_marker[InN_layer,'L'], InN_layer, sep='    ')) )+xlab('')+ylab('')+
  theme(legend.title  = element_blank(),
        axis.text.y = element_text(size = 8,face ='italic'),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

MGE_CGE_layer <-  rev(c("ADARB2" ,"LAMP5", 'VIP','NR2F2',"LHX6", 'SST','PVALB' ,'SOX6','SATB1' ))

InN_mge_order <- plyr::mapvalues(paste('InN',c(6,5,8,1,4,7,3,10,2,9),sep=''),o.ids,n.ids)

h.combined_InN <- subset(h.combined,cluster=='InN')
h.combined_InN@active.ident <- factor(h.combined_InN@active.ident,levels=InN_mge_order)

in_df <- as.data.frame(h.combined_InN[MGE_CGE_layer, ][['RNA']]@data)

in_df$gene <- rownames(in_df)
in_df <- melt(in_df, id = "gene")

in_df$cluster1 <- meta_data[ as.character(in_df$variable),'cluster1' ]

in_df <- in_df %>% group_by(gene, cluster1) %>% dplyr::summarise(meanExp = mean(value)) %>% ungroup
in_df$cluster1 <- factor(in_df$cluster1, InN_mge_order)
in_df$gene <- factor(in_df$gene, levels =  MGE_CGE_layer)

myplot <- ggplot(data = in_df, aes(x=(cluster1), y=(gene), fill=meanExp)) + 
  geom_tile()+ scale_fill_gradientn(colours = myPalette(100))+ coord_equal()+ theme_bw()+
  xlab('')+ylab('')+theme(axis.text.y = element_text(face='italic',size = 8),
                          axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 8))

#############
## Figure 3
#############

# UMAP 
p6_seurat <- subset(h.combined, Period =='P6')
p1 <- DimPlot(p6_seurat,label=T,label.size = 2.7)
p2 <- DimPlot(p6_seurat,label = F,group.by = 'Area') 
myplot <- plot_grid(p1,p2,rel_widths = c(3.5,2.7))
ggsave(filename ='/home1/songlt/Documents/ScDatabase/Figures_corrected/Papers/Fig3_p6_UMAP.png', plot = myplot,width =10 , height = 4.25,units = 'in',dpi=500)

period <- 'P6'
as<- meta_data[meta_data$Period==period,]
as <- as.data.frame(table(as[,c('cluster','Area')]))
as <- subset(as, Freq!=0)

meta_data_s <- as.data.frame(table(subset(meta_data,Period==period)[,c('Area' ,'cluster' )] ))
meta_data_s <- subset(meta_data_s, Freq!=0)
colnames(meta_data_s) <- c('Brain region', 'Celltype','Percentage')

cls <- c("NPC", "ExN" , "InN" , "OPC" , "Olig", "Astro",  "Micro","Endo",  "Perc", "Gran" ,  "Purk"  )

myPlot <- ggplot(meta_data_s, aes(x = `Brain region`, y = Percentage)) + 
  geom_bar(position = "fill",stat = "identity", aes(fill = Celltype)) +
  theme_light()+
  theme(legend.position='top' ,axis.text.x = element_text(vjust = 0.5,size=10.5),title  = element_text(face ='plain'),)+
  
  ggtitle(paste(period,' (',period_infos[period,'age'], ')' ,sep=''))+
  xlab('Brain region')+
  scale_fill_manual(values = c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1'
                               ,'plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue'
  )[which(cls%in%sort(unique(meta_data_s$Celltype)))])

meta_data_s1 <- subset(meta_data_s, `Brain region` %in% c('FC','OC','PC','TC','IG','MDL','Pons'))
colnames(meta_data_s1)[2] <- 'Cell type'
myPlot <- ggplot(meta_data_s1, aes(x = `Brain region`, y = Percentage)) + 
  geom_bar(position = "fill",stat = "identity",width = 0.65) +
  aes(fill =`Cell type` )+
  theme( axis.text.x = element_text(vjust = 0.5,size=10.5),title  = element_text(face ='plain'),legend.key.size = unit(0.3,'cm')) +
  theme_bw()  + 
  xlab('Brain region')+
  scale_fill_manual(values = c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1'
                               ,'plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue'
  )[which(cls%in%sort(unique(meta_data_s$Celltype)))])+   coord_flip()

dev.print(pdf, file='~/Documents/ScDatabase/Figures_corrected/Papers_0522/Fig3_P6_Area_comp2.pdf')

ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/Fig3_P6_Area_comp2.png', plot=myPlot, width = 3, height = 4,units = 'in', dpi=600 )
ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/Fig3_P6_Area_comp3.png', plot=myPlot, dpi=600 )



FC <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'FC','Percentage'], labels='',cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
OC <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'OC','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
PC <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'PC','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
TC <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'TC','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
IG <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'IG','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
MDL <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'MDL','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))
pons <- base2grob(~pie(meta_data_s1[meta_data_s1$`Brain region` == 'Pons','Percentage'],labels = '', cex=0.7, col=c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1','plum2' ,'salmon1' ,'paleturquoise3' ,'paleturquoise4','blue')))

ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_fc_pie.pdf', plot =FC, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_oc_pie.pdf', plot =OC, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_pc_pie.pdf', plot =PC, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_tc_pie.pdf', plot =TC, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_ig_pie.pdf', plot =IG, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_mdl_pie.pdf', plot =MDL, width = 3, height = 3,units = 'in')
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig3_pons_pie.pdf', plot =pons, width = 3, height = 3,units = 'in')

ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_fc_pie.png', plot =FC, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_oc_pie.png', plot =OC, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_pc_pie.png', plot =PC, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_tc_pie.png', plot =TC, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_ig_pie.png', plot =IG, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_mdl_pie.png', plot =MDL, width = 3, height = 3,units = 'in',dpi=600)
ggsave(filename = '~/Documents/ScDatabase/Figures_corrected/Papers/fig3_pons_pie.png', plot =pons, width = 3, height = 3,units = 'in',dpi=600)



CEll_Types_Periods <- function(period){
 
  as<- meta_d[meta_d$Period==period,]
  
  as <- as %>% group_by(Area) %>%
    dplyr::mutate(Proportion = Freq/sum(Freq)) %>% as.data.frame()
  colnames(as)[3] <- 'Subtype'
  

  as1 <- subset(as,Area%in%c('FC','OC','PC','TC','IG','MDL','Pons'))
  
  myPlot <- ggplot(as1, aes(x = Subtype, y = Area,size= Proportion)) +
    geom_point(aes(size= Proportion),show.legend=T,color= 'yellowgreen' ) +
    theme_light()+
    theme(axis.text.x = element_text(angle = 90,hjust=1,vjust=0.5))+
    scale_color_manual(values = c('yellowgreen','tan1'))+
    xlab("Cell subtypes") +
    ylab('')
  
}
CEll_Types_Periods(period='P6')

P6_Seurat <- subset( h.combined, Period=='P6' & protocal=='STRT-seq' & cluster1=='Micro')


P6_Seurat@active.ident <- P6_Seurat$Area

p6_area_marker <- FindAllMarkers(P6_Seurat)
p6_area_marker_exn1_4 <- p6_area_marker
p6_M <- subset(p6_area_marker, avg_logFC>0 )
p6_M <- p6_M[order(p6_M$cluster,-p6_M$avg_logFC),]
marker_freq <- as.data.frame(table(p6_M$gene))
marker_freq <- merge(p6_M, marker_freq, by.x='gene',by.y='Var1',all.x=T,sort=F) 

p6_m_uniq <- as.data.frame(subset(marker_freq, pct.1 >0.5 & pct.2 <0.5)) %>% group_by(cluster) %>% top_n(n = 3, wt =avg_logFC)
p6_m_uniq <- p6_m_uniq[order(p6_m_uniq$cluster) , ]

p6_df <- as.data.frame(P6_Seurat[['RNA']]@data[p6_m_uniq$gene,]) 

anno_col_all <- P6_Seurat@meta.data[order(P6_Seurat@meta.data$Area),c('Area',"Sex",'protocal')]
colnames(anno_col_all) <- c('Region','Gender', 'Platform')
anno_col_all1 <- data.frame(P6_Seurat@meta.data[order(P6_Seurat@meta.data$Area),c('Area')])

rownames(anno_col_all1) <- rownames(anno_col_all)
colnames(anno_col_all1) <- 'Region'


myplot1 <- pheatmap(p6_df[,rownames(anno_col_all)],cluster_rows=F,cluster_cols=F,show_colnames=F,
                    show_rownames=T,fontsize_row=7,annotation_col=anno_col_all      ,
                    legend_breaks=seq(0,5,1),
                    breaks=c(seq(0,1.7,by=0.17),seq(1.71,2.1,0.007),seq(2.11, 3 ,by=0.036),seq(3.01, 4.1, by=0.18) ))



### p6 基因空间表达
Period_geneE <- function(period, gene){
  

  d <- as.matrix(h.combined[['RNA']]@data[gene,])
  
  d <- as.matrix(h.combined[['RNA']]@data[gene,])
  dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
  dd$Gene <- gene
  dd$cell <- as.character(dd$cell)
  dd$cluster <- meta_data[dd$cell,'cluster']
  
  dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
  dd<- dd[!is.na(dd$cluster),]
  dd <- dd[dd$Period==period,]
  
  dd <- subset(dd, Area%in%c('FC','OC','PC','TC','IG','MDL','Pons'))
  p2 <- ggplot(dd,aes(x=Area, y=value, fill=Area, color=Area,alpha=0.4)) +
    geom_violin( size=0.5) + 
    geom_boxplot(width=0.1, outlier.size = 0.1)+
    facet_wrap(~cluster) +
    theme_light()+
    theme(legend.position='' , axis.text.x = element_text(angle = 90,  vjust = 0.5))+
    ylab('Expression Value') +
    xlab('Brain region') +
    ggtitle(paste( c(  period, ' (', gene, ')'), collapse = '',sep=''))
   ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers_0522/Fig3_P6_CD68.pdf',  plot=p2, width = 5, height = 4.5,units = 'in' )
  
  ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/Fig3_P6_CD68.png',  plot=p2, width = 6, height = 4.5,units = 'in', dpi=600 )
  
}


Period_geneE(gene='CD68',period='P6')



CEll_Types_Area <- function(area){
  
  as<- meta_d[meta_d$Area==area,]
  
  as <- as %>% group_by(Period) %>%
    dplyr::mutate(Proportion = Freq/sum(Freq)) %>% as.data.frame()
  
  myPlot <- ggplot(as, aes(x = (cluster1), y = Period)) +
    geom_point(aes(size= Proportion ),color='yellowgreen') +
    theme_light()+
    theme(legend.position = 'top', axis.text.x = element_text(angle = 40,hjust=1))+
    scale_color_manual(values = c('yellowgreen','tan1'))+
    ylab("Period") +
    xlab("Cell Subtype") +
    ggtitle(area)
  
  ggplot(as, aes(x = (Period), y = Proportion, group=cluster1,color=cluster1)) +
    geom_line() +
    theme_light()+
    theme(legend.position = 'top', axis.text.x = element_text(angle = 40,hjust=1))+
    scale_color_manual(values = c('yellowgreen','tan1'))+
    ylab("Period") +
    xlab("Cell Subtype") +
    ggtitle(area)
  
  meta_d <- as.data.frame(table(meta_data[,c('Area','Period','cluster')]))
  meta_d <- meta_d[meta_d$Freq!=0,]
  area <- 'PFC'
  as<- meta_d[meta_d$Area==area,]
  as <- subset(as, !cluster%in%c('Gran','Purk'))
  
  as <- as %>% group_by(Period) %>%
    dplyr::mutate(Proportion = Freq/sum(Freq)) %>% as.data.frame()
  
  myPlot <-  ggplot(as, aes(x = Period, y = Proportion, group=cluster,color=cluster,fill=cluster)) +
    geom_smooth(method = 'loess',span=0.9, se=F)+
    theme_bw()+     theme(legend.title = element_blank())+
    scale_color_manual(values = c('yellow3','tan1','yellowgreen','peachpuff', 'pink1', 'palevioletred1'
                                  ,'plum2' ,'salmon1' ,'paleturquoise3'))
  
  
 
}
CEll_Types_Area('PFC')

### dfc gene expression
Area_geneE <- function(Region, gene){
  
  p1 <- FeaturePlot(object =subset( h.combined, Area==Region), features =  gene )
  
  d <- as.matrix(h.combined[['RNA']]@data[gene,])
  
  dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
  dd$Gene <- gene
  dd$cell <- as.character(dd$cell)
  dd$cluster <- meta_data[dd$cell,'cluster1']
  dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
  dd<- dd[!is.na(dd$cluster),]
  dd <- dd[dd$Area==Region,]
  
  dd$Period <- factor(dd$Period,levels=c(paste('P',1:15,sep='')))
  
  p2 <- ggplot(dd,aes(x=Period, y=value, color=cluster,group=cluster,fill=cluster)) +
    geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15)+
    xlab('Period') +
    ylab('Expression Value')+
    ggtitle(paste(Region," (", gene , ")" , collapse = '',sep=''))+
    theme_light()+
    theme(legend.key.size = unit(0.5, "cm")) 
  
  ggsave(filename=paste('~/Documents/ScDatabase/Figures_corrected/Papers_0522/', 'fig3',Region, '_', gene, '.pdf', sep=''), plot=p2)
  

}

for(Region in c('DFC')) {
  
   for(gene in c('VIP','SYP')){ 
      Area_geneE(gene=gene,Region=Region)
  }
}


#############
## fIGURE 4
#############

## enrichment  fisher test 

Matrix::rowSums(h.combined@assays$RNA@data[1:8000,]==0) -> rs1
Matrix::rowSums(h.combined@assays$RNA@data[8001:16000,]==0) -> rs2
Matrix::rowSums(h.combined@assays$RNA@data[16001:24098,]==0) -> rs3

rs <- c(rs1,rs2,rs3)
#h.combined@assays$integrated  <- NULL
to_inte <- rownames(h.combined)[rs < ncol(h.combined)*0.90]

enrichment_analysis <- function(disease_genes,disease_geneset,marker_df){
  #pv <- data.frame(row.names = unique(marker_df$cluster ) )
  
  for(cl in unique(marker_df$cluster )){                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
    #ge <-marker_df[marker_df$cluster==cl & marker_df$p_val_adj<0.05,'gene']
    ge <-marker_df[marker_df$cluster==cl ,'gene']
    aa <- length(intersect(disease_genes, ge))
    ab <- length(ge) -aa
    ba <- length(disease_genes)
    bb <- length(to_inte) - ba
    
    fr <- fisher.test(matrix(c(aa,ab,ba,bb),nrow = 2,byrow = F),alternative = "greater")
    #pv[cl,disease_geneset]<- fr$p.value
    pv <- rbind(pv,c(disease_geneset, cl, fr$p.value,fr$estimate))
  }
  return(pv)
}

load('~/Documents/ScDatabase/RData/markers_res2_pc30_refine4.RData')

pv <- c()
for (f in dir('./hbst_disease','*txt',full.names = T)){
  
  disease_g <- read.table(f)
  markerfile <- h.combined.markers[h.combined.markers$p_val_adj<0.05,]
  pv <- enrichment_analysis(disease_g$V1,str_split_fixed(f,'\\/|\\.',5)[4], markerfile)
}
pv <- as.data.frame(pv)
colnames(pv) <- c('disease_geneset','cluster','pvalue','OR')

pv$`p-val`  <- as.numeric(as.character(pv$pvalue))
pv$OR <- as.numeric(pv$OR)
pv$cluster <- factor(pv$cluster,levels = levels(meta_data$cluster1))
pv$disease_geneset <- factor(pv$disease_geneset, levels = rev(c('PD','AD','ALS','HD','ASD','BP','MS','Height')))
pv$`-log10(p-val)` <- -log10(pv$`p-val`)
pv1 <- subset(pv, OR!=0)
pv2 <- subset(pv, pvalue<0.05)
colnames(pv)[4] <- 'Odds ratio'
pv$`p-val < 0.05` <- 'N'
pv$`p-val < 0.05`[pv$pvalue<0.05] <- 'Y'
p1 <- ggplot(pv, aes(cluster, disease_geneset)) +
  geom_point()+
  scale_color_distiller(palette = "Spectral")+
  aes(color=`-log10(p-val)` )+
  aes(size=`Odds ratio`)+
  aes(shape=`p-val < 0.05`)+
  theme_bw()+
  theme( legend.key.size=unit(0.4,'cm' ) ,axis.text.x = element_text(angle = 90,vjust = 0.5,hjust=1,size = 7) )+
  xlab('Cell subtype')+
  ylab('Brain disorder')

ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers_0522/Fig4_disorder_enrichment.pdf', plot=p1, width = 9,height = 4)

ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/Fig4_disorder_enrichment.png', plot=p1,dpi = 500,width = 9,height = 4)

asd <- read.table('~/Documents/ScDatabase/processed_data/disease_gene_sets20190906/asd_comb.txt')
h.combined.markers[h.combined.markers$gene%in%asd$V1, ]->aa

set.seed(1234)

samp <- sample(colnames(h.combined),20000)
h.combined_asd <- h.combined[asd$V1,samp]
#h.combined_asd <- h.combined[unique(cell_type_info$C_marker),samp]
asd_e <- as.data.frame(h.combined_asd[['RNA']]@data) 

anno_col_all <- h.combined_asd@meta.data[order(h.combined_asd@meta.data$cluster1),c('cluster1','cluster')]

colnames(anno_col_all) <- c('subtypes','cell types')
anno_col_all$`cell types` <- factor( anno_col_all$`cell types`, levels= unique(anno_col_all$`cell types`))
#asd_g <- c('DSCAM','NRXN1','ANK2','AKAP9','KMT2A','CHD2','ASH1L','TRIO','MYT1L','GRIN2B','GABRB3','SCN2A','SYNGAP1','CHD8')
asd_g <- c('DSCAM','NRXN1','TRIO','SETD5','CTTNBP2','ANK2','AKAP9','KMT2A','CHD2','ASH1L','MYT1L','GRIN2B','GABRB3','SCN2A','SYNGAP1','CHD8')

asd_od <- setdiff(asd$V1,asd_g)
myplot <- pheatmap(asd_e[c(asd_g,asd_od),rownames(anno_col_all)],cluster_rows=F,cluster_cols=F,show_colnames=F,
                   show_rownames=T,fontsize_row=7,annotation_col=anno_col_all,legend_breaks=seq(0,5.6,1),
                   breaks=c(seq(0,2.2,by=0.1),seq(2.21,2.8,0.011),seq(2.81,4.5,by=0.078),seq(4.51,5.1, 0.5 )))

myplot <- pheatmap(asd_e[c(asd_g,asd_od),rownames(anno_col_all)],cluster_rows=F,cluster_cols=F,show_colnames=F,
                   show_rownames=F,fontsize_row=7,annotation_col=anno_col_all,legend_breaks=seq(0,5.6,1),
                   breaks=c(seq(0,4.5,by=0.0478),seq(4.51,5.1, 0.5 )))
ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig4_asd_heatmap3.pdf', plot=myplot, width = 9,height =7)
ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig4_asd_heatmap3_2.pdf', plot=myplot, width = 9,height =15)

ggsave(filename='~/Documents/ScDatabase/Figures_corrected/Papers/fig4_asd_heatmap3.png', plot=myplot,dpi = 500, width = 9,height =7)


## expression of genes associated with brain disorders

Area_geneE <- function(Region, gene){
  
  #d <- as.matrix(h.combined[['RNA']]@data[gene,])
  
  d <- as.matrix(h.combined[['RNA']]@data[gene,])
  
  dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
  dd$Gene <- gene
  dd$cell <- as.character(dd$cell)
  dd$cluster <- meta_data[dd$cell,'cluster1']
  dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
  dd<- dd[!is.na(dd$cluster),]
  dd <- dd[dd$Area==Region,]
  
  dd$Period <- factor(dd$Period,levels=c(paste('P',1:15,sep='')))
  #dd <- subset(dd,cluster%in%c('ExN10', 'ExN9','OPC1','OPC2','OPC3'))
  dd$cluster <- as.character(dd$cluster)
  dd$cluster[!dd$cluster%in%c('ExN9','ExN10', 'OPC1','OPC2','OPC3')] <- 'Others'
  colnames(dd)[4] <- 'Subtype'
  
  P2 <- ggplot(dd,aes(x=Period, y=value, color=Subtype,group=Subtype,fill=Subtype)) +
    geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15)+
    xlab('Period') +
    ylab('Expression Value')+
    ggtitle(paste(Region," (", gene, ")", collapse = '',sep=''))+
    theme_light()+
    theme(legend.key.size = unit(0.5, "cm")) 
  
  
  ggsave(filename=paste('~/Documents/ScDatabase/Figures_corrected/Papers_0522/', 'fig4',Region, '_', gene, '.pdf', sep=''), plot=P2, width =5, height =3   )
  

}

for(Region in c('DFC','PFC')) {
  for(gene in c('DSCAM','NRXN1')){ 
    
    Area_geneE(gene=gene,Region=Region)
  }
}


gene <- 'DSCAM'
#Region <- 'PFC'
Region <- 'DFC'

d <- as.matrix(h.combined[['RNA']]@data[gene,])

dd <- reshape2::melt(d, id = row.names) %>% dplyr::rename(cell = Var1, Gene = Var2)
dd$Gene <- gene
dd$cell <- as.character(dd$cell)
dd$cluster <- meta_data[dd$cell,'cluster1']
dd[,c('Area','Period')] <- meta_data[dd$cell,c('Area','Period')]
dd<- dd[!is.na(dd$cluster),]
dd <- dd[dd$Area==Region,]
dd$cluster <- as.character(dd$cluster)
dd$cluster[!dd$cluster%in%c('ExN9','ExN10', 'OPC1','OPC2','OPC3')] <- 'Others'
tgc <- summarySE(dd, measurevar="value", groupvars=c("cluster","Period"))
pd <- position_dodge(0.1)
# ggplot(tgc, aes(x=Period, y=value, colour=cluster,group=cluster)) + 
#   geom_errorbar(aes(ymin=value-ci, ymax=value+ci), width=.1, position=pd) +
#   geom_line(position=pd) +
#   geom_point(position=pd)

colnames(dd)[4] <- 'subtype'
ggplot(dd,aes(x=Period, y=value, color=subtype,group=subtype,fill=subtype)) +
  geom_smooth(span=1,method = c("auto", "lm", "glm", "gam", "loess")[5],se=T,alpha = 0.15)+
  xlab('Period') +
  ylab('Expression Value')+
  ggtitle(paste(Region," (", gene, ")", collapse = '',sep=''))+
  theme_light()+
  theme(legend.key.size = unit(0.5, "cm")) 

dev.print(pdf, file='~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig4_pfc_dscam.pdf')
dev.print(pdf, file='~/Documents/ScDatabase/Figures_corrected/Papers_0522/fig4_dfc_dscam.pdf')


#############
## Figure S1
#############

a <- as.data.frame(table(meta_data[,c('Area','Period','dataset')]))
b <- as.data.frame(table(meta_data[,c('Area','Period')]))
a <- a[a$Freq!=0,]


p1 <- ggplotly(ggplot(a, aes(x = Period, y = Area, fill = dataset,color=dataset)) +
                 geom_point(aes(size=Freq), shape = 21,alpha=0.65,  position=position_jitter(h=0,w=0.05)) +
                 geom_text(aes(label=''))+theme_bw()+ xlab('Developmental period')+ylab('Brain region')+labs(alpha='',size='Dataset')) 


for(i in 1:nrow(a)){
  a[i,'s'] <- paste(a$dataset[i],' [',a$Freq[i],']',collapse = '',sep='')
}

aa <- dcast(a,Area ~ Period, fun.aggregate=function(x) paste(x, collapse=", ") )
aaa <- melt(aa,id.vars = 'Area')

aaa$value <- gsub(',','\n',aaa$value)
aaa$variable <- factor(aaa$variable,levels = paste('P',1:15,sep=''))

aaa$a <- paste(aaa$Area,aaa$variable,sep='_')
b$a <- paste(b$Area,b$Period,sep='_')
aaab <- merge(aaa,b,by='a')
aaab <- aaab[aaab$Freq!=0,]
colnames(aaab)[2] <- c('Region') 
p1 <- ggplot(aaab,aes(x=Period,y=Region,label = value)) + 
  theme_bw()+
  geom_point(aes(size=Freq,alpha=0.01,color=Region))+
  geom_text(size=3)+theme()+
  xlab('Period')+ylab('Brain region')



#################
##   Figure S3
#################

unique(meta_data$priCluster)[grep('IN|In',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Ex|EX|ex',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('RG|NasN|NEP|IPC|Prog',unique(meta_data$priCluster))]  # hOMTN hNbM 中脑
unique(meta_data$priCluster)[grep('Oli',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Mic',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Peric',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Ast',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Per',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('Gran',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('End',unique(meta_data$priCluster))]
unique(meta_data$priCluster)[grep('OPC',unique(meta_data$priCluster))]


meta_data[,c('cluster1','priCluster')]->bb
bb <- bb[!is.na(bb$priCluster)&bb$priCluster!=''& !grepl('^h',bb$priCluster) & bb$priCluster!='Unk', ]
bb$cluster1 <- as.character(bb$cluster1)
bb$priCluster <- as.character(bb$priCluster)
bb$cluster <- gsub('[0-9]','', bb$cluster1)
bb[grep('IN|In|div',bb$priCluster),'cluster_o'] <- 'InN'
bb[grep('Ex|EX|ex|EN',bb$priCluster),'cluster_o'] <- 'ExN'
bb[grep('RG|NasN|NEP|IPC|Prog',bb$priCluster),'cluster_o'] <- 'NPC'
bb[grep('Oli',bb$priCluster),'cluster_o'] <- 'Olig'
bb[grep('Mic',bb$priCluster),'cluster_o'] <- 'Micro'
bb[grep('Ast',bb$priCluster),'cluster_o'] <- 'Astro'
bb[grep('Per',bb$priCluster),'cluster_o'] <- 'Perc'
bb[grep('End',bb$priCluster),'cluster_o'] <- 'Endo'
bb[grep('OPC',bb$priCluster),'cluster_o'] <- 'OPC'
bb[grep('Gran',bb$priCluster),'cluster_o'] <- 'ExN'
bb[grep('Purk',bb$priCluster),'cluster_o'] <- 'InN'
bb[grep('Gran',bb$cluster),'cluster'] <- 'ExN'
bb[grep('Purk',bb$cluster),'cluster'] <- 'InN'

bb$concordancy <- c(bb$cluster==bb$cluster_o)
bb<- bb[!is.na(bb$concordancy),]
table(bb$concordancy)[2]/nrow(bb)

bb <- dcast(as.data.frame(table(bb[,c('cluster','concordancy')])),cluster~concordancy)
bb$Concordance <- bb$`TRUE`/(  bb$`FALSE`+bb$`TRUE`  )
bb$percent <- paste(round(bb$Concordance*100,1),'%',sep='')
#table(bb[,c('concordancy')])

#barplot(bb$Concordance,names=bb$cluster ,col="#69b3a2",labels)
ggplot(data = bb, aes(x=cluster,y=Concordance))+
  geom_bar(stat = "identity",fill="#69b3a2")+
  geom_text(data=bb,aes(x=cluster,y=Concordance,label=percent),vjust=0) +
  xlab('Cell types')

dev.print(pdf, file='/home1/songlt/Documents/ScDatabase/Figures/Papers/cell_type_concoedance.pdf')

#################
##   Figure S3 
#################
#h.combined$cluster <- gsub('\\d','',h.combined$cluster1)
cls <- c("NPC", "ExN" , "InN" , "OPC" , "Olig", "Astro",  "Micro","Endo",  "Perc", "Gran" ,  "Purk"  )
h.combined$cluster <- factor(h.combined$cluster,level=cls)

p1 <- DimPlot(h.combined,group.by = 'cluster1',label = T)+NoLegend()

p2 <- DimPlot(h.combined,group.by = 'dataset',label = F)

p3 <- DimPlot(h.combined,group.by = 'Period',label = F)

p4 <- DimPlot(h.combined,group.by = 'Area',label = F)

#################
##   Figure S5
#################

p4 <- DimPlot(h.combined,group.by = 'cluster1',label = F,split.by = 'dataset')+NoLegend()


#################
##   Figure S6 
#################

lake_marker <- read.table('~/Documents/ScDatabase/processed_data/lake_2015_marker.txt',sep='\t',header = T)


lake_in <-  lake_marker[lake_marker$cluster %in%c(paste('In', 1:8,sep='')) ,'Gene']

h.combined_InN <- subset(h.combined,cluster=='InN')

myplot <- DotPlot(h.combined_InN, features = rev(lake_in)) + RotatedAxis() +
  theme(legend.position='top' , axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.75)))


lake_ex <-  unique(lake_marker[lake_marker$cluster %in%c(paste('Ex', 1:8,sep='')) ,'Gene'])
h.combined_ExN <- subset(h.combined,cluster=='ExN' )

myplot <- DotPlot(h.combined_ExN, features = rev(lake_ex)) + RotatedAxis() +  
  theme(legend.position='' , axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.3)))




#################
##   Figure S7 Inn EXN markers
#################

load('~/Documents/ScDatabase/RData/markers_res2_pc30_refine3.RData')

h.combined.markers <- subset(h.combined.markers, p_val_adj<0.05 )

marker_freq <- as.data.frame(table(h.combined.markers$gene))
marker_freq <- merge(h.combined.markers, marker_freq, by.x='gene',by.y='Var1',all.x=T) 
#marker_freq <- as.data.frame(marker_freq %>% group_by(cluster) %>% top_n(n = 3, wt = avg_logFC ))
final_markers <- as.data.frame(subset(marker_freq, Freq<=2) %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC ))

final_markers <- final_markers[grepl('ExN|InN', final_markers$cluster),] %>% arrange(cluster)

myplot <- plot_vln(subset(h.combined,cluster%in%c('ExN','InN')), gene=final_markers$gene)




#################
##   Figure S8 
#################
num_p <- as.data.frame(table(meta_data_c$Period))
rownames(num_p) <- num_p$Var1
dc <- reshape2::dcast(meta_data_c,Period~cluster1)

for(ap in unique(dc$Period)){
  dc[dc$Period==ap,c(-1)] <- dc[dc$Period==ap,c(-1)]/(num_p[ap,'Freq']/num_c)
}

dc_mp <- melt(dc,id=c('Period'))
dc_mp <- dc_mp[dc_mp$value!=0,]
p_c<- as.data.frame(cbind(levels(dc_mp$Period), c(brewer.pal(n = 12, name = c("Set3")),brewer.pal(n = 8, name = c("Set2")))[1:15]))
rownames(p_c) <- p_c$V1
dc_mp <- merge(dc_mp, p_c, by.x='Period',by.y='V1') 
dc_mp_exn <- dc_mp[grep('ExN',dc_mp$variable),]
p1 <- ggplot(dc_mp_exn, aes(x = variable, y = value)) + 
  geom_bar(position = "fill",stat = "identity", aes(fill = Period)) +
  theme_light()+
  xlab('Excitatory cell subtypes')+
  ylab('Proportion')



#################
##   Figure S9-10 Astro - Olig
#################

plot_vln <- function(t, gene) {
  d <- as.matrix(t[['RNA']]@data[intersect(gene, rownames(t[['RNA']]@data)), ])
  dd <- melt(d, id = row.names)
  dd <- dd %>% dplyr::rename(gene = Var1, cell = Var2)
  dd$tree.ident <- t$cluster1[dd$cell]
  #str(dd$tree.ident)
  dd$gene <- factor(dd$gene, levels = intersect(gene, rownames(t[['RNA']]@data)))
  ggplot(dd, aes(tree.ident, value, fill = tree.ident)) + geom_violin(scale = "width", trim = T, alpha = 0.8, adjust = 1) + 
    facet_wrap(~gene, scales = "free_y", ncol = 1, strip.position = "right") +
    theme(strip.background = element_blank(), strip.placement = "outside", 
          strip.text.y = element_text(colour = "red", angle = 360, size = 10),  panel.grid = element_blank(), axis.title.y = element_blank(),
          panel.border = element_blank()) + 
    theme(strip.text.y = element_text(colour = "red", angle = 360, size = 10,face ='italic'),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = rel(0.9)),
          axis.text.y=element_blank(),
          axis.ticks.y =element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          legend.position="none",
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())+
    scale_y_continuous(limits=c(0,5))
}

markers <- c('AQP4',final_markers[final_markers$cluster=='Astro1', 'gene'][1:3], 
             final_markers[final_markers$cluster=='Astro2', 'gene'][c(2,4,5)],
             final_markers[final_markers$cluster=='Astro3', 'gene'][c(1,2,3)], 
             final_markers[final_markers$cluster=='Astro4', 'gene'][c(1, 2,3)])

myplot <- plot_vln(t=h.combined, gene=markers)



markers <- c('MBP',final_markers[final_markers$cluster=='Olig1', 'gene'][1:2],
             final_markers[final_markers$cluster=='Olig2', 'gene'][3:4],
             final_markers[final_markers$cluster=='Olig3', 'gene'][2:3],
             final_markers[final_markers$cluster=='Olig4', 'gene'][c(3,4)])

myplot <- plot_vln(t=h.combined, gene=markers)


## 画Astro Olig subtype的饼图
meta_data_c <- meta_data
meta_data_c$Area_p <- paste(meta_data_c$Area, meta_data_c$Period,sep='_')

num_a <- as.data.frame(table(meta_data_c$Area))
rownames(num_a) <- num_a$Var1
dc <- reshape2::dcast(meta_data_c,Area~cluster1)

for(ap in unique(dc$Area)){
  dc[dc$Area==ap,c(-1,-2)] <- dc[dc$Area==ap,c(-1,-2)]/(num_a[ap,'Freq']/num_c)
}

dc_ma <- melt(dc,id=c('Area'))
dc_ma <- dc_ma[dc_ma$value!=0,]

for(cl in c(paste('Astro',1:4,sep=''), paste('Olig',1:4,sep=''))){
  dc_mas <- dc_ma[dc_ma$variable==cl,]
  label_a <- paste(dc_mas$Area,' (', round(100*dc_mas$value/sum(dc_mas$value),1),'%',')',sep='')
  a1 <-  base2grob(~pie(dc_mas$value,labels = label_a, cex=0.7, col=c(brewer.pal(n = 12, name = c("Set3")),brewer.pal(n = 12, name = c("Set2")))))
  ggsave(filename = paste('~/Documents/ScDatabase/Figures_corrected/Papers/',cl,'_area.png',sep=''), plot =a1, width = 5, height = 5,units = 'in',dpi=100)
}

#############################
##   Figure S11 ##############
#############################

cl_p <- as.data.frame(xtabs(~cluster1+Period, meta_data))
cl_p <- subset(cl_p, Freq!=0)
#cl_p$Period <- as.character(cl_p$Period)
cl_p_n <- as.data.frame(table(cl_p$Period))
cl_p_n <- subset(cl_p_n, Freq!=0)
myplot <- ggplot(data=cl_p, aes(x=Period,y=Freq))+
  geom_bar(position = "fill",stat = "identity", aes(fill = cluster1))+
  geom_text(data=cl_p_n, aes(Var1, y=1.021, label=Freq), 
            vjust=1, fontface=4, col="black", inherit.aes=FALSE)+
  theme(legend.position = 'top',legend.key.size = unit(0.27, "cm"),legend.text=element_text(size=6.5),legend.title = element_blank())+
  guides(fill=guide_legend(nrow=3))+
  xlab('')


cl_a <- as.data.frame(xtabs(~cluster1+Area, meta_data))
cl_a <- subset(cl_a, Freq!=0)
#cl_a$Period <- as.character(cl_a$Period)
cl_a_n <- as.data.frame(table(cl_a$Area))
cl_a_n <- subset(cl_a_n, Freq!=0)
myplot <- ggplot(data=cl_a, aes(x=Area,y=Freq))+
  geom_bar(position = "fill",stat = "identity", aes(fill = cluster1))+
  geom_text(data=cl_a_n, aes(Var1, y=1.021, label=Freq), 
            vjust=1, fontface=4, col="black", inherit.aes=FALSE)+
  theme(legend.position = 'no',legend.key.size = unit(0.27, "cm"),legend.text=element_text(size=6.5),legend.title = element_blank(),
        axis.text.x = element_text(angle = 30,vjust = 0.7))+
  xlab('')

#############################
##   Figure S12 ##############
#############################

Go_Kegg <- function(cl,marker_df){
  
  #ge <- marker_df[marker_df$cluster==cl & marker_df$p_val_adj<0.05,'gene']
  ge <- marker_df[marker_df$cluster==cl , 'gene']
  
  ge <- bitr(ge, fromType="SYMBOL", toType=c("ENTREZID","ENSEMBL"), OrgDb="org.Hs.eg.db")
  ge <- ge$ENTREZID[!duplicated(ge$ENTREZID)]
  go <- enrichGO(ge, OrgDb = org.Hs.eg.db, ont='BP',pAdjustMethod = 'BH',pvalueCutoff = 0.05, 
                 qvalueCutoff = 0.2,keyType = 'ENTREZID')
  
  gor <- go@result
  gor$cluster <- cl
  
  kegg <- enrichKEGG(ge, organism = "hsa",pvalueCutoff=0.01)
  barplot(kegg,showCategory=20,drop=T)
  dotplot(kegg, showCategory=30)
  
  keggr <- kegg@result
  keggr$cluster <- cl
  
  return(list(go=gor,kegg=keggr))
}

go_df <- c()
kegg_df <- c()

p6_area_marker <- p6_M[p6_M$p_val_adj<0.05, ]
p6_area_marker <- p6_area_marker[order( p6_area_marker[,7], -p6_area_marker[,2] ) ,]

for(cl in p6_area_marker$cluster%>%unique()){
  Go_Kegg_r <-  Go_Kegg(cl,p6_area_marker)
  go_df <- rbind(go_df,Go_Kegg_r$go)
  kegg_df <- rbind(kegg_df,Go_Kegg_r$kegg)
}

go_df <- go_df[!duplicated(go_df$ID), ]
kegg_df <- kegg_df[!duplicated(kegg_df$ID), ]

go_t5 <- go_df %>% group_by(cluster) %>% top_n(5, -pvalue)%>%as.data.frame()
kegg_t5 <- kegg_df %>% group_by(cluster) %>% top_n(5, -pvalue)%>%as.data.frame()
go_kegg_t5 <- rbind(kegg_t5, go_t5)
go_kegg_t5 <- go_kegg_t5[ order(go_kegg_t5$cluster,go_kegg_t5$pvalue),]
go_kegg_t5$Description <- factor(go_kegg_t5$Description, levels=go_kegg_t5$Description)


myplot <- ggplot(data=go_kegg_t5,aes(x=Description,y=-log10(pvalue),fill = cluster),fill = cluster)+
  geom_bar(stat = 'identity')+
  coord_flip()+
  theme_bw()+
  xlab('')+
  scale_fill_discrete(name = "")

