if (T) {
  dir.create("scripts")
  dir.create("results")
  dir.create("files")
  dir.create("figures")
  dir.create("origin_datas/GEO",recursive = T)
  dir.create("origin_datas/TCGA")
}

library(Seurat)
library(harmony)
library(tidyverse)
library(sctransform)
library(cols4all)
library(data.table)
library(stringr)
library(ComplexHeatmap)
library(scales)
library(ggsci)

dir.create("01_landscape")

assays <- dir("origin_datas/GEO/GSE139829_RAW/")
dir <- paste0("origin_datas/GEO/GSE139829_RAW/", assays)
assays
dir

samples_name = assays
samples_name

scRNAlist <- list()
for(i in 1:length(dir)){
  counts <- Read10X(data.dir=dir[i])
  scRNAlist[[i]] <- CreateSeuratObject(counts, project=samples_name[i], min.cells=3, min.features=200)
  scRNAlist[[i]] <- RenameCells(scRNAlist[[i]], add.cell.id=samples_name[i])   
  if(T){    
    scRNAlist[[i]][["percent.mito"]] <- PercentageFeatureSet(scRNAlist[[i]], pattern="^MT-") 
  }
}


names(scRNAlist) <- samples_name

scRNA <- merge(scRNAlist[[1]], scRNAlist[2:length(scRNAlist)])



mycolor <- c4a('poly.alphabet2',26)#c4a("brewer.paired", 12)
show_col(mycolor)
colors <- mycolor
#p = VlnPlot(scRNA, features=c("nFeature_RNA","nCount_RNA",'percent.mito'), pt.size=0.0, cols=colors)
VlnPlot_nFeature_RNA_before = VlnPlot(scRNA,features=c("nCount_RNA",'nFeature_RNA','percent.mito'),pt.size = 0)
VlnPlot_nFeature_RNA_before
ggsave("01_landscape/VlnPlot_nFeature_RNA_before.pdf", VlnPlot_nFeature_RNA_before, 
       width=12, height=5)



scRNA = subset(scRNA, subset=nFeature_RNA>200&nFeature_RNA<8000&percent.mito<10)
VlnPlot_nFeature_RNA_after = VlnPlot(scRNA,features=c("nCount_RNA",'nFeature_RNA','percent.mito'),pt.size = 0)
VlnPlot_nFeature_RNA_after
ggsave("01_landscape/VlnPlot_nFeature_RNA_after.pdf", VlnPlot_nFeature_RNA_after, 
       width=12, height=5)
NormalizeData, FidnVariableFeatures, ScaleData 
# s.genes=Seurat::cc.genes.updated.2019$s.genes
# g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
# scRNA <- CellCycleScoring(scRNA, s.features=s.genes, g2m.features=g2m.genes, set.ident=TRUE)
scRNA = SCTransform(scRNA, vars.to.regress=c( "percent.mito"), verbose=FALSE)
scRNA = RunPCA(scRNA, verbose=FALSE)
scRNA = RunHarmony(scRNA, group.by.vars="orig.ident", max.iter.harmony=50, lambda=0.5, assay.use="SCT")
ElbowPlot(scRNA,ndims = 50)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims=1:30, reduction="harmony")
# p = DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=1)+
#   theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)



save(scRNA,file = 'scRNA.Rdata')
# load('scRNA.Rdata')


ElbowPlot(scRNA,ndims = 50)
scRNA <- FindNeighbors(scRNA, dims=1:30, reduction="harmony")
scRNA <- RunUMAP(scRNA, dims = 1:30,reduction="harmony")
p1 = DimPlot(scRNA, reduction="umap", group.by="orig.ident", pt.size=1)+
  theme(legend.position="right", plot.title=element_blank())+scale_color_manual(values=colors)
ggsave("01_landscape/UMAP_Sample.pdf", p1, height=6, width=6)
resolution=0.1################################ 
#c4a_gui()
colors <- sample(c4a('poly.sky24',24),size = 20,replace = F)
mydata <- FindClusters(scRNA, resolution=0.1)
# p2 <- UMAPPlot(mydata, pt.size=1, label=T, cols=colors, label.size=5)+NoLegend()
# ggsave("01_landscape/cluster_umap.pdf",p2,height = 10,width = 10)
markers_cell_type <- FindAllMarkers(mydata,group_by='cell_type',only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)
write.table(markers_cell_type, "01_landscape/All_cluster_type.txt", col.names=T, row.names=T, quote=F, sep="\t")
table(mydata$seurat_clusters)


VlnPlot(mydata, features=c("MLANA","PMEL","DCT","TYR","TYRP1","MITF","ALCAM","KIT","SOX10","MMP12","PHLDA1","APOD","ATP1B1","CDH19","MX2","FILIP1L","AHNAK2","PTGDS","CRABP1","DKK3","BNC2"), pt.size=0,group.by ='seurat_clusters')+NoLegend()+theme(axis.title.x=element_blank())
#0 Melanocyte cells:MLANA,PMEL
#1 Melanocyte cells:MLANA,PMEL
#2:CD8+ T cells:CD8A CD8B CD3D GZMA
#3 Melanocyte cells:MLANA,PMEL
#4: Macrophage /Monocyte cells	:CD14 CD68 LYZ
#5 Melanocyte cells:MLANA,PMEL
#6: Plasma cells:CD27 CD38 MZB1 IGHG1
#7 Melanocyte cells:MLANA,PMEL
#8 Melanocyte cells:MLANA,PMEL
#9: Myeloid cell:PTGDS
#10 Natural killer T (NKT) cells 1 :KPNA2 CDC25B CENPM IDH2
#11 Endothelial cells ：PECAM1	 VWF CD34	
#12 Melanocyte cells:MLANA,PMEL
#13  Cancer stem cell：PROM1 
#14 Photoreceptor cells 2  ：RCVRN

Idents(mydata) <- mydata$seurat_clusters
table(mydata$SCT_snn_res.0.1)
#mydata = subset(mydata, seurat_clusters %in% c(0,1,2,3,4,6,7,8,9))
#mydata@meta.data$seurat_clusters = droplevels(mydata@meta.data$seurat_clusters, exclude=c(5, 8, 9, 11, 12, 13))

cell_label = c("Melanocyte cells", "Melanocyte cells", "CD8+ T cells", "Melanocyte cells", "Macrophage /Monocyte cells",
                "Melanocyte cells","Plasma cells", "Melanocyte cells","Melanocyte cells",'Melanocyte cells',
               'Melanocyte cells',    'Endothelial cells','Melanocyte cells','Cancer stem cells','Melanocyte cells')

names(cell_label) <- levels(mydata)
mydata <- RenameIdents(mydata, cell_label)
mydata[["cell_type"]] = Idents(mydata)
# c4a_gui()
colors = c4a('tableau.20',19)[c(1,3,5,7,9,11,13)]

p2 = UMAPPlot(mydata, pt.size=1, label=T, label.size=6)+NoLegend()+scale_color_manual(values=colors)+theme(text = element_text())  #family = "Arial"
 
library(ggpubr)
p2=ggarrange(AugmentPlot(p2,dpi = 300,width = 12,height = 8), ncol = 1, nrow = 1,
             legend ='right',legend.grob=get_legend(p2),
             font.label = list(size = 6, face = "bold",family ='Times')) 
p2
ggsave("01_landscape/UMAP_cell_type2.pdf", p2, width=10, height=8)
genes = c('MLANA','PMEL' ,#'TYRP1','MITF','ALCAM',#'TYR',
          'CD8A','CD3D', 'CD8B', 'GZMA',
          'CD14', 'CD68', 'LYZ',
           'CD38', 'MZB1', 'IGHG1',#'CD27',
          #'PTGDS',
         #'KPNA2', 'CDC25B', 'CENPM', 'IDH2',
          'PECAM1',	'CD34',	'VWF',
          'PROM1'#'BEX1'
          #'RCVRN'
          )

#0 Melanocyte cells:MLANA,PMEL
#1 Melanocyte cells:MLANA,PMEL
#2: T cells:CD8A CD8B CD3D GZMA
#3 Melanocyte cells:MLANA,PMEL
#4: Macrophage /Monocyte cells	:CD14 CD68 LYZ
#5 Melanocyte cells:MLANA,PMEL
#6: Plasma cells:CD27 CD38 MZB1 IGHG1
#7 Melanocyte cells:MLANA,PMEL
#8 Melanocyte cells:MLANA,PMEL
#9: Myeloid cell:PTGDS
#10 Natural killer T (NKT) cells 1 :KPNA2 CDC25B CENPM IDH2
#11 Endothelial cells ：PECAM1	 VWF CD34	
#12 Melanocyte cells:MLANA,PMEL
#13  Cancer stem cell：PROM1 
#14 Photoreceptor cells 2  ：RCVRN




       
p3 = DotPlot(mydata, features=genes)+coord_flip()+scale_color_gradientn(colors=c("#211C84", "#7A73D1", "white", "#B3D8A8", "#3D8D7A"))+theme_minimal()+
  theme(text=element_text(family="Times"), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), axis.title.x=element_blank(), axis.title.y=element_blank())
p3
ggsave("01_landscape/dotplot_gene_marker.pdf",p3,height = 8,width = 10)



genes2 <-c('MLANA','PMEL' ,#'TYRP1','MITF','ALCAM',#'TYR',
           'CD8A','CD3D', #'CD8B', 'GZMA',
           'CD14', 'CD68', 'LYZ',
           'CD38', 'MZB1', #'IGHG1',#'CD27',
           #'PTGDS',
           #'KPNA2', 'CDC25B',# 'CENPM', 'IDH2',
           'PECAM1',	'CD34',	#'VWF',
           'PROM1'#'BEX1'
           #'RCVRN'
)
p3_2 <- VlnPlot(mydata, features=genes2, pt.size=0, cols=colors,group.by ='cell_type',ncol = 5)+NoLegend()+theme(axis.title.x=element_blank())
p3_2 
ggsave("01_landscape/volin_marker.pdf", p3_2, width=12, height=8)
# p3_3 <- VlnPlot(mydata, features=c("MLANA","PMEL","DCT","TYR","TYRP1","MITF","ALCAM","KIT","SOX10","MMP12","PHLDA1","APOD","ATP1B1","CDH19","MX2","FILIP1L","AHNAK2","PTGDS","CRABP1","DKK3","BNC2"), pt.size=0)+NoLegend()+theme(axis.title.x=element_blank())
# p3_3
# ggsave("01_landscape/dotplot_gene_marker_3.pdf",p3_3,height = 25,width = 25) 

save(mydata,file = 'mydata.Rdata')




mydata[['site']]=ifelse(mydata$orig.ident%in%c('GSM4147091','GSM4147092','GSM4147100'),'metastatic','primary')
Type_label = c("metastatic", "primary")
bar <-  as.data.frame(with(mydata@meta.data, table(site, cell_type)))
bar$site = factor(bar$site, levels=Type_label)
bar = bar %>% group_by(site) %>% mutate(percent=100*Freq/sum(Freq))

p4 <- ggplot(data=bar, aes(x=cell_type, y=percent, fill=site,label = sprintf("%.2f", percent)))+
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_manual(values=c('#6A9C89','#FFA725'))+theme_classic()+
  ggtitle("Percent(%)")+
  geom_text(position = position_dodge(width = 0.9), vjust = -0.5, size = 4)+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_text(angle=30, hjust=1, size=12, face="bold"), axis.text.y=element_text(face="bold", size=12), legend.text=element_text(face="bold", size=12), legend.title=element_blank())
p4 
ggsave("01_landscape/barplot_pair_number.pdf", p4, width=8, height=4)

library(lessR)
mata.data.all <- mydata
celltype.all <-  as.data.frame(mata.data.all@meta.data[,'cell_type'])
colnames(celltype.all) <-'cell_type' 
rownames(celltype.all) <-rownames(mata.data.all@meta.data)
pdf("01_landscape/all_PieChart.pdf",height =8,width = 12)
PieChart( cell_type,data = celltype.all,
         hole =0.5,
         main="Cell Proportion of all",
         main_cex =1.3,fill=colors
         )
dev.off()


