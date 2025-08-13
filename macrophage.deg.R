
library(stringr)
library(tidydr)
library(openxlsx)
library(data.table)
library(reshape2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(clusterProfiler)
library(pheatmap)
library(ComplexHeatmap)
library(GSVA)
library(GSEABase)
library(fgsea)
library(corrplot)
library(colorspace)
library(survival)
library(survminer)
library(maftools)
library(vegan)
library(forcats)
library(ggpubr)
library(ggsci)
library(ggplot2)
library(rstatix)
library(ggstatsplot)
library(ggcor)
library(ggstance)
options(stringsAsFactors = F)
source('/home/pub252/projects/codes/mg_base.R')
mg_violin_1=function(data,xangle=0,ylab='value',xlab='',leg.title='Group',test_method='anova',cmp_test_method='t.test',
                     legend.pos='r',melt=F,jitter=T,ylim=NULL,show_compare=NULL,point_size=NULL,group_col){
  library(ggplot2)
  if(is.null(ylim)){
    
  }
  if(melt){
    data_m=data
    colnames(data_m)=c('Group','value')
  }else{
    data_m=reshape2::melt(data)
    colnames(data_m)=c('Group','value')
  }
  if(!is.null(ylim)){
    data_m$value[data_m$value>ylim[2]]<-NA
  }
  data_m=data_m[which(!is.na(data_m[,1])),]
  if(xangle==0){
    tx=element_text(colour="black",family="Times")
  }else{
    tx=element_text(angle=xangle,hjust = 1,colour="black",family="Times")
  }
  
  pos='right'
  if(is.null(legend.pos)){
    pos='none'
  }else if(legend.pos=='tr'){
    pos=c(1,1)
  }else if(legend.pos=='br'){
    pos=c(1,0)
  }else if(legend.pos=='tl'){
    pos=c(0,1)
  }else if(legend.pos=='bl'){
    pos=c(0,0)
  }else if(legend.pos=='t'){
    pos='top'
  }else if(legend.pos=='r'){
    pos='right'
  }else if(legend.pos=='b'){
    pos='bottom'
  }
  uni.group=unique(data_m[,1])
  ct=length(uni.group)
  
  p1<-ggplot(data_m,aes(x=Group,y=value))+geom_violin(alpha=0.7)
  # if(ct<=4){
  #   p1=p1+ggsci::scale_fill_lancet()
  # }else if(ct<=10){
  #   p1=p1+ggsci::scale_fill_npg(name=leg.title)
  # }else if(ct<=20){
  #   p1=p1+ggsci::scale_fill_d3(palette = "category20",name=leg.title)
  # }else if(ct<=30){
  #   cbPalette=c(ggsci::pal_npg("nrc", alpha = 0.6)(10),ggsci::pal_d3("category20", alpha = 0.6)(20))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }else if(ct<=38){
  #   cbPalette=c(ggsci::pal_lancet()(10)
  #               ,ggsci::pal_npg("nrc", alpha = 0.6)(10)
  #               ,ggsci::pal_d3("category20", alpha = 0.6)(20)
  #               ,ggsci::pal_nejm("default", alpha = 0.6)(8))
  #   p1=p1+scale_fill_manual(values=cbPalette[1:ct])
  # }
  p1=p1+scale_fill_manual(values=group_col)
  # if(jitter){
  #   if(is.null(point_size)){
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2)
  #   }else{
  #     p1<-p1+geom_jitter(alpha=0.3,col='black',show.legend=FALSE,width = 0.2,size=point_size)
  #   }
  # }
  # 
  p1=p1+theme_bw()+geom_boxplot(width=0.2,aes(fill=Group),outlier.shape = NA)
  p1=p1+theme(axis.text.x=tx, 
              axis.text.y=element_text(family="Times",face="plain"), 
              axis.title.y=element_text(family="Times",face="plain"), 
              legend.text=element_text(face="plain", family="Times", colour="black" 
              ),
              legend.title=element_text(face="plain", family="Times", colour="black"
              ),
              legend.justification=pos, legend.position=pos
              ,legend.background = element_rect(fill = NA, colour = NA)
  )+ylab(ylab)+xlab(xlab)
  til=''
  if(test_method=='anova'){
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=t.test(x1,x2)$p.value 
      til=paste0('t-tests p=',signif(pv,2))
    }else{
      fit <- aov(value~Group, data = data_m)
      pv=summary(fit)[[1]][5][[1]]
      fv=summary(fit)[[1]][4][[1]]
      til=paste0('ANOVA tests p=',signif(pv,2))
    }
  }else{
    if(length(unique(data_m[,1]))<3){
      x1=data_m[,2][which(data_m[,1]==unique(data_m[,1])[1])]
      x2=data_m[,2][which(data_m[,1]==unique(data_m[,1])[2])]
      pv=wilcox.test(x1,x2)$p.value 
      til=paste0('wilcox.tests p=',signif(pv,2))
    }else{
      fit=kruskal.test(value~Group, data = data_m)
      pv=fit$p.value
      til=paste0('Kruskal-Wallis test p=',signif(pv,2))
    }
  }
  p1=p1+ggtitle(til) 
  if(!is.null(ylim)){
    p1=p1+ylim(ylim)
  }
  if(is.null(show_compare)){
    if(length(uni.group)>5){
      show_compare=F
    }else{
      show_compare=T
    }
  }
  if(show_compare){
    comps=list()
    for(i in 1:(length(uni.group)-1)){
      for(j in (i+1):length(uni.group)){
        comps=c(comps,list(c(uni.group[i],uni.group[j])))
      }
    }
    p1=p1+ggpubr::stat_compare_means(comparisons = comps,method = cmp_test_method,label= "p.signif")
  }
  return(p1)
}
my_mutiboxplot=function(dat,group,group_cols=ggsci::pal_aaas()(10),
                        #test_method=c('t.test','wilcox.test','anova','kruskal.test')[4],
                        bw=T,xlab='',ylab='score',title='',size=10,angle = 45, hjust = 1,
                        legend.position='top',fill='group',notch=F){
  # dat=tcga.est[tcga.subtype.cli$Samples,]
  # group=tcga.subtype.cli$Cluster
  dat.bind=cbind(dat,Cluster=group)
  dat.bind=crbind2DataFrame(dat.bind)
  dat.melt=melt(dat.bind)
  #data=data[which(!is.na(data[,1])),]
  colnames(dat.melt)=c('Group','type','value')
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=dat.melt %>%
    ggplot(aes(x=type, y=value,fill=Group)) +
    geom_boxplot(notch = notch) +  
    scale_fill_manual(values =group_cols)+   
    ggpubr::stat_compare_means(aes(group=Group), label = "p.signif", method = test_method)+
    labs(x="", y = ylab, fill =fill,title =title) +
    #theme_light()+
    theme_bw()+
    #theme_classic()
    theme(legend.position = legend.position,                
          plot.title = element_text(hjust = 0.5),text = element_text(family = 'Times',size = size),
          axis.text.x = element_text(angle = angle, hjust = hjust),
          panel.grid = element_blank()) 
  return(p)
}
my_violin=function(dat,group,group_cols=ggsci::pal_aaas()(10),#test_method='kruskal.test',
                   fill= "Group",label=c("p.format",'p.signif')[1],
                   xlab='',ylab='',title='',x.size=10,y.size=10,legend.position='top'){
  
  # dat = tcga.b.cell$GeneSet,
  # group = tcga.subtype$Cluster
  data=data.frame(Group=group,value=dat)
  data=crbind2DataFrame(data)
  data=melt(data)
  data=data[which(!is.na(data[,1])),]
  if(length(names(table(group)))>2){
    test_method='kruskal.test'
  }else{
    test_method='wilcox.test'
  }
  p=ggplot(data,aes(x=Group, y=value,fill=Group)) +
    geom_violin(trim = F)+  
    scale_fill_manual(values = group_cols)+
    geom_boxplot(width=0.2,position=position_dodge(0.8),outlier.colour = NA,fill="white")+
    theme_classic(base_size = 20)+labs(x=xlab,y=ylab,title = title)+
    ggpubr::stat_compare_means(aes(group=Group), label = label, method =test_method)+
    theme(legend.position = legend.position,axis.text = element_text(color = 'black'),
          title = element_text(size = 12),
          axis.text.x = element_text(size = 12),axis.title.x = element_text(size = 15),
          axis.text.y = element_text(size = 12),axis.title.y = element_text(size = 12))
  return(p)
}
custom_theme <- function() {
  theme_survminer() %+replace%
    theme(text = element_text(family = 'Times'),panel.grid = element_blank())
}
getGeneFC=function(gene.exp,group,ulab=ulab,dlab=dlab){
  degs_C1_C3=mg_limma_DEG(gene.exp, 
                          group,
                          ulab=ulab,
                          dlab=dlab)
  
  degs_C1_C3_sig<-degs_C1_C3$DEG[which(degs_C1_C3$DEG$adj.P.Val <= 1),]
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(stringr)
  degs_C1_C3_sig_gene<-rownames(degs_C1_C3_sig)
  
  degs_gene_entz=bitr(degs_C1_C3_sig_gene,fromType="SYMBOL",toType="ENTREZID",OrgDb="org.Hs.eg.db")
  degs_gene_entz <- dplyr::distinct(degs_gene_entz,SYMBOL,.keep_all=TRUE)
  
  gene_df <- data.frame(logFC=degs_C1_C3_sig$logFC,
                        SYMBOL = rownames(degs_C1_C3_sig))
  gene_df <- merge(gene_df,degs_gene_entz,by="SYMBOL")
  head(gene_df)
  
  geneList<-gene_df$logFC
  names(geneList)=gene_df$ENTREZID 
  head(geneList)
  
  geneList=sort(geneList,decreasing = T)
  head(geneList) 
  return(geneList)
}
get_riskscore.lasso<-function(dat,os,os.time,labels=c('A','B')){
  library(glmnet)
  # set.seed(2024)
  fit1=glmnet(as.matrix(dat)
              ,cbind(time=os.time,
                     status=os)
              ,family="cox"
              ,nlambda=100
              , alpha=1) 
  
  cv.fit<-cv.glmnet(as.matrix(dat)
                    ,cbind(time=os.time,
                           status=os)
                    ,family="cox"
                    ,nfolds = 10
                    ,nlambda=100
                    , alpha=1)
  sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
  #print(cv.fit$lambda.min)
  #length(names(sig.coef))
  #10
  mg_plot_lasso <- function(fit,cv_fit,lambda=NULL,show_text=T,figLabels=c('A','B')){
    if(is.null(lambda)){
      lmda=cv_fit$lambda.min
    }else{
      lmda=lambda
    }
    fit.coef=fit$beta[(apply(fit$beta,1,function(x){
      return(sum(x!=0))
    })>0),]
    
    fit.coef=as.matrix(fit.coef)
    colnames(fit.coef)=fit$lambda
    #fit$lambda==cv_fit$lambda
    library(ggplot2)
    dat=data.table::melt(t(as.matrix(fit.coef)))
    dat_z=dat[which(dat$value==0),]
    dat=dat[which(dat$value!=0),]
    dat.sv=rbind()
    for (u in unique(dat_z[,2])) {
      t.z=dat_z[which(dat_z[,2]==u),1]
      t.zx=max(t.z)
      dat.sv=rbind(dat.sv,c(t.zx,u,0))
      t.zn=min(t.z)
      if(t.zx!=t.zn){
        dat.sv=rbind(dat.sv,c(t.zn,u,0))
      }
    }
    colnames(dat.sv)=colnames(dat_z)
    #dat_z=dat_z[dat_z[,2]%in%names(which(fit.coef[,which(fit$lambda==lmda)]!=0)),]
    dat=crbind2DataFrame(rbind(dat,dat.sv))
    mn=min(-log(dat$Var1))
    mx=max(-log(dat$Var1))
    if(show_text){
      mx=(mx-mn)*0.1+mx
    }
    p=ggplot(dat, aes(x=-log(Var1), y=value,colour=Var2))+geom_line()+theme_bw()+theme(legend.position = "none")
    p=p+coord_cartesian(xlim=c(mn, mx))+xlab('-ln(lambda)')+ylab('Coefficients')
    if(show_text){
      fl=fit.coef[which(fit.coef[,which(fit$lambda==lmda)]!=0),ncol(fit.coef)]
      for_label=data.frame(Var1=rep(min(dat$Var1),length(fl)),Var2=names(fl),value=fl)
      p=p+ggrepel::geom_label_repel(
        aes(label = Var2,color=Var2),
        data = for_label,hjust = 0
      )
    }
    p=p+geom_vline(aes(xintercept=-log(lmda)), colour="#BB0000", linetype="dashed")
    p=p+annotate('text',x=-log(lmda),y=min(dat[,3]),label=paste0('lambda=',round(lmda,4)))
    tgc=data.frame(lambda=cv_fit$lambda,cvm=cv_fit$cvm,cvup=cv_fit$cvup,cvlo=cv_fit$cvlo,cvsd=cv_fit$cvsd
                   ,col=ifelse(cv_fit$lambda>=cv_fit$lambda.min&cv_fit$lambda<=cv_fit$lambda.1se,ifelse(cv_fit$lambda==lmda,'A','C'),'B'))
    p1=ggplot(tgc, aes(x=log(lambda), y=cvm)) + xlab('ln(lambda)')+ ylab('Partial Likelihood Deviance')+
      geom_errorbar(aes(ymin=cvm-cvsd, ymax=cvm+cvsd)) +
      geom_point(aes(colour=col))
    p1=p1+theme_bw()+theme(legend.position = "none")
    gal=ggpubr::ggarrange(p,p1, ncol = 2, nrow = 1
                          #,align = "hv"
                          ,labels = figLabels)
    return(gal)
  }
  lasso.pdf <- mg_plot_lasso(fit1,
                             cv.fit,
                             show_text=T,
                             figLabels=labels)
  return(list(lasso.gene=names(sig.coef),lambda.min=cv.fit$lambda.min,plot=lasso.pdf))
}
my_riskplot=function(cli_dat,cols=c("red","blue"),xlab='Samples',
                     a.ylab="Risk score",b.labs="Survival time(year)",cutoff=0,labs=c('A','B')){
  #cli_dat=tcga.risktype.cli
  cli.dat.order=cli_dat[order(cli_dat$Riskscore),c('OS.time','Status','Riskscore','Risktype')]
  fp_dat=data.frame(Samples=1:nrow(cli_dat),cli.dat.order)
  p1=ggplot(fp_dat,aes(x=Samples,y=Riskscore))+geom_point(aes(color=Risktype))+
    scale_colour_manual(values =cols)+
    theme_bw()+labs(x=xlab,y=a.ylab)+
    geom_hline(yintercept=cutoff,colour="black", linetype="dotted",size=0.8)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p1
  p2=ggplot(fp_dat,aes(x=Samples,y=OS.time))+geom_point(aes(col=Status))+theme_bw()+
    scale_colour_manual(values =cols)+
    labs(x=xlab,y=b.labs)+
    geom_vline(xintercept=sum(fp_dat$Risktype=="Low"),colour="black", linetype="dotted",size=0.8)
  p2
  p.merge=mg_merge_plot(p1,p2,nrow=2,ncol=1,labels = labs)
  return(p.merge)
}

library(openxlsx)
library(Seurat)
library(dplyr)
library(ggplot2)
library(magrittr)
library(gtools)
library(stringr)
library(Matrix)
library(tidyverse)
library(patchwork)
library(data.table)
library(RColorBrewer)
library(ggpubr)
library(dplyr)
library(tidydr)
library('ggsci')
#install.packages("SeuratTheme") 
library(SeuratTheme) 
library(gridExtra)
library(cowplot)
#devtools::install_github('junjunlab/scRNAtoolVis')
library(scRNAtoolVis)
tcga.data=read.delim('origin_datas/TCGA/Merge_RNA_seq_FPKM (1).txt',row.names = 1,check.names = F)
range(tcga.data)
tcga.data[1:4,1:4]
#table(substr(colnames(tcga.data),14,15))
dim(tcga.data)
#80
# tcga.data$gene <- rownames(tcga.data) 
# tcga.data$gene  <- sub("\\..*", "", tcga.data$gene )
# tcga.data <- tcga.data %>%
#   group_by(gene) %>%
#   summarise(across(everything(), mean, na.rm = TRUE)) %>%
#   column_to_rownames("gene")  
tcga.data <- exp_ensg2symbol(tcga.data)
genecode=read.delim('data/GeneTag.genecode.v32.txt',sep='\t',header = T)
table(genecode$TYPE)
mrna_genecode=genecode[which(genecode$TYPE=='protein_coding'),]
head(mrna_genecode)

range(tcga.data)
tcga.exp=log2(tcga.data[rownames(tcga.data) %in% mrna_genecode$SYMBOL,]+1)
range(tcga.exp)
#colnames(tcga.exp) <- substr(colnames(tcga.exp), 1, nchar(colnames(tcga.exp)) - 1)

tcga.pancancer.cli=read.xlsx('data/TCGA_pancancer_cli_PMID_29625055.xlsx')
head(tcga.pancancer.cli)
tcga.cli=tcga.pancancer.cli[which(tcga.pancancer.cli$type=='UVM'),]
head(tcga.cli)
tcga.cli=data.frame(Samples=paste0(tcga.cli$bcr_patient_barcode,'-01'),
                    Age=tcga.cli$age_at_initial_pathologic_diagnosis,
                    Gender=tcga.cli$gender,
                    AJCC_stage=tcga.cli$ajcc_pathologic_tumor_stage,
                    tcga.cli[,c('OS','OS.time','DSS','DSS.time','DFI','DFI.time','PFI','PFI.time')])
rownames(tcga.cli)=tcga.cli$Samples
head(tcga.cli)
tcga.cli$PFI.time
#tcga.cli=tcga.cli %>% drop_na(PFI.time)
tcga.cli=tcga.cli[tcga.cli$OS.time>0,]
dim(tcga.cli)
head(tcga.cli)
fivenum(tcga.cli$Age)
tcga.cli$Age1=ifelse(tcga.cli$Age>60,'>61.5','<=61.5')
table(tcga.cli$Gender)
table(tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage[tcga.cli$AJCC_stage=='[Not Available]']=NA
tcga.cli$AJCC_stage=gsub('[ABC]','',tcga.cli$AJCC_stage)
tcga.cli$AJCC_stage=gsub('Stage ','',tcga.cli$AJCC_stage)
tcga.exp <- tcga.exp[,intersect(tcga.cli$Samples,colnames(tcga.exp))]
tcga.cli <- tcga.cli[intersect(tcga.cli$Samples,colnames(tcga.exp)),]
tcga.exp <- tcga.exp[,intersect(tcga.cli$Samples,colnames(tcga.exp))]
dim(tcga.exp)
########GSE22138
# library(GEOquery)
# GSE22138 <- getGEO(filename = "origin_datas/GEO/GSE22138_series_matrix.txt.gz", getGPL = F)
load('origin_datas/GEO/GSE22138.RData')
GSE22138.pheno=pData(GSE22138)
head(GSE22138.pheno)
table(GSE22138.pheno$`metastasis:ch1`)
GSE22138.cli=data.frame(Samples=GSE22138.pheno$geo_accession,
                        event=GSE22138.pheno$`metastasis:ch1`,
                        age=GSE22138.pheno$`age:ch1`,
                        gender=GSE22138.pheno$`gender:ch1`,
                        time=as.numeric(GSE22138.pheno$`months to endpoint:ch1`)/12*365)
rownames(GSE22138.cli)=GSE22138.pheno$Samples
GSE22138.cli$event=ifelse(GSE22138.cli$event=='no',0,1)
head(GSE22138.cli)


GSE22138.exp=exprs(GSE22138)
GSE22138.exp[1:5,1:5]
GSE22138.exp=exp_probe2symbol_v2(GSE22138.exp,GPL = 'GPL570')
range(GSE22138.exp)
dim(GSE22138.exp)




dir.create("02_DEG_Plasma")
all_cell_deg <- read.table("01_landscape/All_cluster_type.txt",sep = "\t",header = T,row.names = 1)

cell_type_mapping <- list(
'0'="Melanocyte cells", 
'1'="Melanocyte cells", 
'2'="CD8+ T cells", 
'3'="Melanocyte cells", 
'4'="Macrophage /Monocyte cells",
'5'="Melanocyte cells",
'6'="Plasma cells", 
'7'="Melanocyte cells",
'8'="Melanocyte cells",
'9'='Melanocyte cells',
'10'='Melanocyte cells', 
'11'='Endothelial cells',
'12'='Melanocyte cells',
'13'='Cancer stem cells',
'14'='Melanocyte cells'
)
library(dplyr)

all_cell_deg$cluster <- sapply(all_cell_deg$cluster, function(x) {
  cell_type_mapping[[as.character(x)]]
})
colnames(all_cell_deg)[2] <- 'avg_log2FC'
# all_cell_deg <- all_cell_deg %>%
#   rename(
#     'avg_log2FC' = avg_logFC
#     #cluster1 = 6,
#    # 'cluster1'= 'cluster',
#    # 'cell_type'='cluster'
#   )
#save(all_cell_deg,file = "02_DEG_Plasma/all_cell_deg.Rdata")
all_cell_figure <- jjVolcano(diffData = all_cell_deg,
                             #tile.col = colors[1:2],
                             pvalue.cutoff = 0.05,
                             log2FC.cutoff = 0.25,
                             topGeneN = 0)
ggsave("02_DEG_Plasma/all_cell_figure.pdf",all_cell_figure,height = 8,width = 18)
 write.table(all_cell_deg,file = "01_landscape/All_Cell_type_DEG.txt",sep = "\t",quote = F,col.names = T,row.names = T)


table(all_cell_deg$cluster)
deg.Plasma <- all_cell_deg[all_cell_deg$cluster=="Plasma cells",'gene']
wgcna.gene.enrichment=mg_clusterProfiler(genes = deg.Plasma )
head(wgcna.gene.enrichment$GO_BP)
wgcna_enrich_res=rbind(wgcna.gene.enrichment$KEGG@result,wgcna.gene.enrichment$GO_BP@result,wgcna.gene.enrichment$GO_CC@result,wgcna.gene.enrichment$GO_MF@result)
wgcna_enrich_res$DB <- rep(c('KEGG','GO_BP','GO_CC','GO_MF'),c(nrow(wgcna.gene.enrichment$KEGG@result),nrow(wgcna.gene.enrichment$GO_BP@result),nrow(wgcna.gene.enrichment$GO_CC@result),nrow(wgcna.gene.enrichment$GO_MF@result)))


write.xlsx(
  list(GO_BP=wgcna.gene.enrichment$GO_BP,
       KEGG=wgcna.gene.enrichment$KEGG,GO_CC=wgcna.gene.enrichment$GO_CC,GO_MF=wgcna.gene.enrichment$GO_MF),'02_DEG_Plasma/wgcna.gene.enrichment.xlsx',
  overwrite = T)
colnames(wgcna_enrich_res)
wgcna_enrich_res <- wgcna_enrich_res %>%
  group_by(DB) %>%
  filter(p.adjust > 0.05) %>%   
  arrange(DB)   
wgcna_enrich_res$Description=factor(wgcna_enrich_res$Description
                                    ,levels=wgcna_enrich_res$Description[order(wgcna_enrich_res$DB,decreasing = T)], ordered=TRUE)
head(wgcna_enrich_res)
fig2b<- ggplot(data = wgcna_enrich_res, aes(x = Description, y = -log10(p.adjust), fill = DB, group = Description)) +
  geom_bar(stat="identity", position="dodge", colour=aes(DB))+
  scale_fill_manual(values =pal_simpsons()(9)[5:7])+
  theme_classic()+
  theme(text = element_text(family = 'Times',size = 16),axis.text.x = element_text(vjust = 0,angle = 90))
fig2b
ggsave("02_DEG_Plasma/wgcna_go.pdf",fig2b,height =10,width = 12 )

##MCells.ssgsea####
library(survival)
library(survminer)
deg.Plasma <- as.character(deg.Plasma)
MCells.score <- t(ssGSEAScore_by_genes(gene.exp = tcga.exp,genes =deg.Plasma ))

#03.WGCNA
dir.create('03.WGCNA')
library(WGCNA)
allowWGCNAThreads(nThreads = 36)
enableWGCNAThreads(nThreads = 36)

my_mad <- function(x){mad(x,na.rm = TRUE)} 
wgcna_exp=t(tcga.exp)
m.mad <- apply(wgcna_exp,2,my_mad)
dim(tcga.exp)
# 

tpm_T2 <- wgcna_exp[,which(m.mad >max( quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01))]


range(tpm_T2)

pdf('03.WGCNA/1.pdf',width = 8,height = 8)
tpm_T2.power=mg_wgcna_get_power(tpm_T2)
dev.off()


tpm_T2.power$cutPower
tpm_T2.module=mg_WGCNA_getModule(tpm_T2,
                                 power = tpm_T2.power$cutPower,
                                 deepSplit=2,
                                 mergeCutHeight=0.25,
                                 minModuleSize=60)

table(tpm_T2.module$Modules[,2])
length(table(tpm_T2.module$Modules[,2]))

pdf('03.WGCNA/2.pdf',height = 5,width = 12)
plotDendroAndColors(tpm_T2.module$Tree, tpm_T2.module$Modules,
                    c("Dynamic Module",'Merged Module'),
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()
MODULE.gene.num <- as.data.frame(table(tpm_T2.module$Modules[,2]))
write.table(MODULE.gene.num,file = "03.WGCNA/MODULE.gene.num.txt",sep = "\t",quote = F,row.names =F) 
writeMatrix(tpm_T2.module$Modules,outpath = '03.WGCNA/tcga.wgcna.module.genes.txt')
pdf('03.WGCNA/3.pdf',height = 6,width = 6)
mg_barplot_point(labels = names(table(tpm_T2.module$Modules[,2]))
                 ,values = as.numeric(table(tpm_T2.module$Modules[,2]))
                 ,point_sizes = 2
                 ,point_cols = names(table(tpm_T2.module$Modules[,2]))
                 ,xlab = 'Number of Genes',legend.pos = NULL)
dev.off()


# Calculate eigengenes
MEs = tpm_T2.module$MEs
# Calculate dissimilarity of module eigengenes
MEDiss = 1-cor(MEs);
# Risktype module eigengenes
METree = hclust(as.dist(MEDiss), method = "average");
# Plot the result
pdf('03.WGCNA/4.pdf',height = 6,width = 12,onefile = T)
plot(METree, main = "Risktypeing of module eigengenes",xlab = "", sub = "")
dev.off()

tcga_cli_use <-data.frame(MCells.score=MCells.score[tcga.cli$Samples,])
head(tcga_cli_use)
tcga_cli_use.part=tcga_cli_use
str(tcga_cli_use.part)

tcga_cli_use.part=sapply(tcga_cli_use.part, function(x)as.numeric(as.factor(x)))
spms=tcga_cli_use.part

MEs_col<-tpm_T2.module$MEs
dim(MEs_col)

modTraitCor = cor(MEs_col[,rownames(MEDiss)[METree$order]]
                  , spms
                  ,use = 'pairwise.complete.obs')
modTraitP = corPvalueStudent(modTraitCor, dim(spms)[1])

textMatrix = paste(signif(modTraitCor, 2), " (", format(modTraitP,scientific =TRUE,digits = 3), ")", sep = "")
dim(textMatrix) = dim(modTraitCor)
dim(textMatrix)

rownames(modTraitCor)=gsub("ME","",rownames(modTraitCor))
rownames(textMatrix)=gsub("ME","",rownames(textMatrix))
colnames(modTraitCor)

pdf('03.WGCNA/5.pdf',width =4,height = 10)
labeledHeatmap(Matrix = data.frame(modTraitCor),
               xLabels = colnames(modTraitCor),
               yLabels = rownames(modTraitCor),
               cex.lab = 1,
               ySymbols = colnames(t(modTraitCor)), colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = data.frame(textMatrix),
               setStdMargins = FALSE,xLabelsAngle = 0,
               cex.text = 0.8, zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


geneModuleMembership <- signedKME(tpm_T2
                                  , data.frame(tpm_T2.module$MEs)
                                  , outputColumnName = "")
head(geneModuleMembership)
MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership)
                                           , nrow(tpm_T2.module$MEs)))

geneTraitSignificance <- as.data.frame(cor(tpm_T2
                                           , spms
                                           , use = 'pairwise.complete.obs'))

GSPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance)
                                           , nrow(spms)))

modNames<-colnames(geneModuleMembership)
modNames

module = "brown"
column = match(module, modNames)
column
moduleGenes <- (tpm_T2.module$Modules[,'mergedColors']==module)
tcga.wgcna.genes=names(which(moduleGenes))
length(tcga.wgcna.genes)
#1063

pdf('03.WGCNA/6.pdf',height = 8,width = 8)
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, "MCells.score"]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for MCells.score",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2
                   , col = module,lwd=2)
dev.off()

dir.create('04.model')
tcga.wgcna.genes <- deg.Plasma
#length(tcga.wgcna.genes)
#1063

tcga.wgcna.genes.cox=cox_batch(tcga.exp[tcga.wgcna.genes,tcga.cli$Samples],
                             time = tcga.cli$PFI.time,
                             event = tcga.cli$PFI)
table(tcga.wgcna.genes.cox$p.value<0.05)
tcga.wgcna.genes.cox=na.omit(tcga.wgcna.genes.cox)
pre.genes=rownames(tcga.wgcna.genes.cox[tcga.wgcna.genes.cox$p.value<0.05,])

bioForest(rt = tcga.wgcna.genes.cox[pre.genes,],col = c('#EB5B00','#89AC46'))



tcga_model_data <- cbind(tcga.cli[, c("OS.time", "OS")],
                         t(tcga.exp[pre.genes, tcga.cli$Samples]))
colnames(tcga_model_data) <- gsub('-', '_', colnames(tcga_model_data))

##LASSO########
library(glmnet)
set.seed(2023)
fit1=glmnet(as.matrix(tcga_model_data[,-c(1,2)])
            ,cbind(time=tcga_model_data$OS.time,
                   status=tcga_model_data$OS)
            ,family="cox"
            ,nlambda=100
            , alpha=1) 

cv.fit<-cv.glmnet(as.matrix( tcga_model_data[,-c(1,2)])
                  ,cbind(time=tcga_model_data$OS.time,
                         status=tcga_model_data$OS)
                  ,family="cox"
                  ,nfolds = 10
                  ,nlambda=100
                  , alpha=1)

sig.coef <- coefficients(cv.fit,s=cv.fit$lambda.min)[which(coefficients(cv.fit,s=cv.fit$lambda.min)[,1]!=0),1]
print(cv.fit$lambda.min)
length(names(sig.coef))
names(sig.coef)

pdf('04.model/LASSO.pdf',height = 3.5,width = 6,onefile = F)
par(mfrow=c(1,2))
plot(fit1)
plot(cv.fit)
dev.off()


fmla <- as.formula(paste0("Surv(OS.time,OS) ~"
                          ,paste0(names(sig.coef),collapse = '+')))
cox <- coxph(fmla, data =as.data.frame(tcga_model_data))
 cox=step(cox)
lan <- coef(cox)
lan
paste0(round(lan, 3), '*', names(lan),collapse = '+')
"0.732*ISG20+-2.095*PDE4B+0.899*SDF2L1"

gene.coef=data.frame(gene=names(lan),coef=as.numeric(lan))
gene.coef
gene.coef$Type=ifelse(gene.coef$coef>0,'Risk','Protective')
gene.coef.fig=ggplot(gene.coef, aes(x = coef, y = reorder(gene,coef), fill =Type)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = c("#FFA725", "#6A9C89")) +
  labs(x = 'coefficient', y = "") +
  geom_text(aes(label = round(coef,3),hjust =2), data = subset(gene.coef, coef > 0))+ 
  geom_text(aes(label = round(coef,3), hjust = -1), data = subset(gene.coef, coef < 0))+  
  theme_bw()+ theme(text = element_text(family = 'Times',size = 14,face = 'bold'),legend.position = 'top')
gene.coef.fig
ggsave('04.model/gene.coef.fig.pdf',gene.coef.fig,height = 3.5,width = 4)

module.coxforest=ggforest(cox, data = tcga_model_data, 
                          main = "Hazardratio", fontsize =1.0, 
                          noDigits = 2)
module.coxforest
ggsave('04.model/module.coxforest.pdf',module.coxforest,height = 3.5,width = 6)




risktype.col=c('#D98324','#443627')
risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)


# tcga.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
#                             variables = 'Riskscore')
# tcga.point.cutoff <- as.numeric(summary(tcga.point)[1])
# tcga.point.cutoff
# tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.point.cutoff,'High','Low')

tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>median(tcga.risktype.cli$Riskscore),'High','Low')

tcga.roc.OS=ggplotTimeROC(tcga.risktype.cli$OS.time,
                           tcga.risktype.cli$OS,
                           tcga.risktype.cli$Riskscore,mks = c(1,2,3))
tcga.roc.OS
tcga.km.OS=ggsurvplot(fit=survfit(Surv(OS.time/365, OS) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='TCGA-UVM',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Progression Free Interval(OS)',
                       legend.position='top',
                       ggtheme = theme_bw(base_size = 12))
tcga.km.OS=mg_merge_plot(tcga.km.OS$plot,tcga.km.OS$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.OS


# risk.tcga=as.numeric(lan%*%as.matrix(t(tcga_model_data[tcga.cli$Samples,names(lan)])))
# tcga.risktype.cli=data.frame(tcga.cli,Riskscore=risk.tcga)


# tcga.point <- surv_cutpoint(tcga.risktype.cli, time = "OS.time", event = "OS",
#                             variables = 'Riskscore')
# tcga.point.cutoff <- as.numeric(summary(tcga.point)[1])
# tcga.point.cutoff
# tcga.risktype.cli$Risktype=ifelse(tcga.risktype.cli$Riskscore>tcga.point.cutoff,'High','Low')



tcga.roc.PFI=ggplotTimeROC(tcga.risktype.cli$PFI.time,
                           tcga.risktype.cli$PFI,
                           tcga.risktype.cli$Riskscore,mks = c(1,2,3))
tcga.roc.PFI
tcga.km.PFI=ggsurvplot(fit=survfit(Surv(PFI.time/365, PFI) ~ Risktype,
                                   data = tcga.risktype.cli),
                       data=tcga.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='TCGA-UVM',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Progression Free Interval(PFI)',
                       legend.position='top',
                       ggtheme = theme_bw(base_size = 12))
tcga.km.PFI=mg_merge_plot(tcga.km.PFI$plot,tcga.km.PFI$table,nrow=2,heights = c(3,1),align = 'v')
tcga.km.PFI



##GSE22138
GSE22138_model_data <- data.frame(GSE22138.cli[, c("time", "event")],
                                  t(GSE22138.exp[intersect(names(lan),rownames(GSE22138.exp)), GSE22138.cli$Samples]))
colnames(GSE22138_model_data) <- gsub('-', '_', colnames(GSE22138_model_data))

intersect(names(lan),rownames(GSE22138.exp))
risk.GSE22138=as.numeric(lan%*%as.matrix(t(GSE22138_model_data[GSE22138.cli$Samples,names(lan)])))
GSE22138.risktype.cli=data.frame(GSE22138.cli,Riskscore=risk.GSE22138)


GSE22138.risktype.cli$Risktype=ifelse(GSE22138.risktype.cli$Riskscore>median(GSE22138.risktype.cli$Riskscore),'High','Low')
GSE22138.roc=ggplotTimeROC(GSE22138.risktype.cli$time,
                           GSE22138.risktype.cli$event,
                           GSE22138.risktype.cli$Riskscore,mks = c(1,2,3))
GSE22138.roc
GSE22138.km=ggsurvplot(fit=survfit(Surv(time/365, event) ~ Risktype,
                                   data = GSE22138.risktype.cli),
                       data=GSE22138.risktype.cli,
                       conf.int = T,pval = T,risk.table = T, 
                       fun = "pct",size = 1,surv.median.line = 'hv',
                       title='GSE22138',legend.title='Risktype',
                       legend.labs = c('High','Low'),
                       legend.position='top',
                       linetype = c("solid", "dashed","strata")[1],
                       palette = risktype.col,ylab='Progression Free Interval(PFI)',
                       ggtheme = theme_bw(base_size = 12))
GSE22138.km=mg_merge_plot(GSE22138.km$plot,GSE22138.km$table,nrow=2,heights = c(3,1),align = 'v')
GSE22138.km





dir.create('05_Risktype.immune')

# ##estimate####
tcga_estimate <- immu_estimate(exp = tcga.exp)
head(tcga_estimate)
my_mutiboxplot(dat = tcga_estimate[tcga.risktype.cli$Samples,],group = tcga.risktype.cli$Risktype,group_cols = risktype.col)

tcga_estimate <- tcga_estimate[,1:3]
tcga_estimate <- cbind(tcga_estimate[tcga.risktype.cli$Samples,],
                       tcga.risktype.cli)

colnames(tcga.risktype.cli)
tcga_StromalScore_cor <- cor_point(x = tcga_estimate$StromalScore,
                                   y = tcga_estimate$Riskscore,
                                   xlab = 'StromalScore',
                                   ylab = 'RiskScore',top_col='#3D8D7A',right_col='#A4B465'
)
tcga_StromalScore_cor

tcga_ImmuneScore_cor <- cor_point(x = tcga_estimate$ImmuneScore,
                                  y = tcga_estimate$Riskscore,
                                  xlab = 'ImmuneScore',
                                  ylab = 'RiskScore',
                                  top_col='#8EB486',right_col='#997C70')
tcga_ImmuneScore_cor

tcga_ESTIMATEScore_cor <- cor_point(x = tcga_estimate$ESTIMATEScore,
                                    y = tcga_estimate$Riskscore,
                                    xlab = 'ESTIMATEScore',
                                    ylab = 'RiskScore',
                                    top_col='#8EB486',right_col='#997C70')
tcga_ESTIMATEScore_cor



p6a<- cowplot::plot_grid(tcga_StromalScore_cor,
                         tcga_ImmuneScore_cor,
                         tcga_ESTIMATEScore_cor,
                         ncol = 3,labels = c("",'',''))
p6a
ggsave("05_Risktype.immune/p5a.pdf",p6a,height = 6,width = 12)
####CIBERSORT#####
tcga.TIMER=read.csv('data/infiltration_estimation_for_tcga.csv',
                    check.names = F,row.names = 1)
tcga.TIMER[1:5,1:5]
tcga.TIMER=tcga.TIMER[tcga.risktype.cli$Samples,]
dim(tcga.TIMER)
table(str_split_fixed(colnames(tcga.TIMER),'_',2)[,2])


tme.df2=tcga.TIMER[tcga.risktype.cli$Samples,str_split_fixed(colnames(tcga.TIMER),'_',2)[,2]=='CIBERSORT']
colnames(tme.df2)=gsub('_CIBERSORT','',colnames(tme.df2))
tme.df2$Risktype=tcga.risktype.cli$Risktype
tme.df2=melt(tme.df2)
head(tme.df2)
p5b<- ggplot(tme.df2,aes(x=variable,y=value,fill=Risktype))+
  geom_boxplot()+stat_compare_means(aes(group=Risktype), label = "p.signif", method = 'wilcox.test')+
  scale_fill_manual(values =risktype.col)+
  xlab('')+ylab('Fraction')+ggtitle('CIBERSORT')+
  theme_bw()+theme(text = element_text(family = 'Times',size = 12),legend.position = 'top',
                   axis.text.x = element_text(color = "black", size = 12,angle = 30,hjust = 1),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())
ggsave('05_Risktype.immune/Fig5c.pdf',height = 5,width = 12)

p5ab<- mg_merge_plot(p6a,p5b,common.legend =T,widths = c(1,1.),labels =c ("A","B"),legend = "top",nrow=2)
savePDF(p5ab,file = "05_Risktype.immune//p5ab.pdf",height = 12,width = 12)
###mcp###
tcga.mcp <- immu_MCPcounter(exp =tcga.exp,isTCGA = T)
saveRDS(tcga.mcp,file ='05_Risktype.immune//taga.mcp.RDS')
mg_PlotMutiBoxplot(data = tcga.mcp[tcga.risktype.cli$Samples,],group =tcga.risktype.cli$Risktype,
                   legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'ssgsea Immune Score')
library(tidyverse)
library(ggcor)
library(vegan)
cr=psych::corr.test(x=tcga.risktype.cli[,'Riskscore'],
                    y=tcga.mcp[tcga.risktype.cli$Samples,]
                    ,method = 'spearman')
df=t(rbind(cr$r,cr$p))
colnames(df)=c('r','p.value')
df=data.frame(Riskscore='Riskscore',MCP_count=rownames(df),df)
df
df <- df %>%
  mutate(lty = cut(r, breaks = c(-1, 0, 1),
                   labels = c("r <= 0", "r > 0")),
         col = cut(p.value, breaks = c(0, 0.01, 0.05, 1),
                   labels = c("< 0.01", "< 0.05", ">= 0.05"),
                   right = FALSE, include.lowest = TRUE))
head(df)
corrmat.color=colorRampPalette(c('#626F47', 'white','#FFCF50'))(100)

p5c<-quickcor(tcga.mcp[tcga.risktype.cli$Samples,], cor.test = TRUE,type = "lower") + #upper
  geom_square(data = get_data(p.value < 0.05, type = "lower")) + 
  anno_link(df, mapping = aes(colour = col,
                              size = abs(r),
                              linetype = lty),diag.label = TRUE) +
  scale_fill_gradient2n(colours = corrmat.color) +
  remove_x_axis()+
  scale_size_area(max_size = 2) +
  scale_linetype_manual(values = c("dotted", "solid")) +
  guides(
    fill = guide_colourbar(title = "corr", order = 3),
    colour = guide_legend(title = "spearman's p", order = 2),
    size = guide_legend(title = "spearman's r", order = 1),
    linetype = "none")
p5c
ggsave('05_Risktype.immune/p5b.pdf',p5c,height = 6,width = 6)




# tcga.exp.icg=immu_ICGs(tcga.exp)
# dim(tcga.exp.icg)
# tcga_icg_plot <- mg_PlotMutiBoxplot(tcga.exp.icg[rownames(tcga.risktype.cli), ]
#                                     , group = tcga.risktype.cli$Risktype
#                                     , legend.pos = 'top'
#                                     , add = 'boxplot'
#                                     , xangle = 60
#                                     , ylab = 'EXP'
#                                     , group_cols = risktype.col
#                                     #, test_method = 'kruskal.test'
# )
# tcga_icg_plot
# ggsave("05_Risktype.immune//tcga_icg_plot.pdf",tcga_icg_plot,height = 6,width = 14)




library(ComplexHeatmap)
library(dplyr)

colnames(tcga.risktype.cli)
tcga.risktype.cli <- arrange(tcga.risktype.cli, Risktype)

tcga_ssgsea_icg <- t(scale(tcga.exp.icg[rownames(tcga.risktype.cli), ]))
pdf('05_Risktype.immune//tcga_ssgsea_icg_plot.pdf', width = 8, height = 8)
tcga_ssgsea_icg_plot <- Heatmap(tcga_ssgsea_icg
                                , name = "Expression"
                                , col = circlize::colorRamp2(c(-2, 0, 2), c( '#003092', 'white','#A31D1D'))
                                , border = T
                                , show_column_names = F
                                , show_column_dend = F
                                , show_row_dend = F
                                , cluster_columns = T
                                , cluster_rows = T
                                , column_split = factor(tcga.risktype.cli$Risktype)
                                , clustering_distance_rows  ='pearson'
                                , clustering_method_rows = 'ward.D2'
                                # , column_km =2
                                , row_names_gp = gpar(fontsize = 10)
                                , top_annotation = HeatmapAnnotation(Risktype = tcga.risktype.cli$Risktype
                                                                     , col=list(Risktype=c('High'=risktype.col[1],
                                                                                           'Low' = risktype.col[2])
                                                                     )
                                                                     , annotation_width = unit(c(1,2), 'cm')
                                                                     , annotation_height = unit(0.2, "cm")
                                                                     , gap = unit(1, 'mm')))
tcga_ssgsea_icg_plot
dev.off()

#6.GSEA###########
dir.create('06_GSEA')
tcga.geneList=getGeneFC(gene.exp=tcga.exp,group=tcga.risktype.cli$Risktype
                        ,ulab='High',dlab = 'Low')
saveRDS(tcga.geneList,file = "06_GSEA/tcga.geneList.RDS")
library(clusterProfiler)
tcga.gsea.KEGG=gseKEGG(
  tcga.geneList,
  organism = "hsa",
  keyType = "kegg",
  pvalueCutoff = 0.05,
  pAdjustMethod = "BH",
  verbose = TRUE,
  use_internal_data = FALSE,
  seed = FALSE,
  by = "fgsea")
saveRDS(tcga.gsea.KEGG,file = "06_GSEA/tcga.gsea.KEGG.RDS")
tcga.gsea.KEGG <- readRDS("06_GSEA/tcga.gsea.KEGG.RDS")

tcga.gsea.KEGG.res=tcga.gsea.KEGG@result
write.csv(tcga.gsea.KEGG.res,'06_GSEA/GSEA_res.csv',row.names = F)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES<0)
table(tcga.gsea.KEGG.res$p.adjust<0.05 & tcga.gsea.KEGG.res$NES>0)
library(dplyr)
ind1=tcga.gsea.KEGG.res %>% slice_max(n =5, order_by = NES)
ind2=tcga.gsea.KEGG.res %>% slice_min(n =5, order_by = NES)  


FIG6a=enrichplot::gseaplot2(tcga.gsea.KEGG,ind1$ID, pvalue_table = F,title ='KEGG enrichment in High group')
FIG6b=enrichplot::gseaplot2(tcga.gsea.KEGG,ind2$ID, pvalue_table = F,title ='KEGG enrichment in Low group')
gsea.p=mg_merge_plot(FIG6a,FIG6b,nrow = 2,labels = c('A','B'))
ggsave('06_GSEA/risktype_GSEA.pdf',gsea.p,height = 15,width = 12)


dir.create("07_drug")
# predictedPtype_Cisplatin <- pRRopheticPredict(as.matrix(tcga_exp)
#                                               , "Cisplatin"
#                                               , selection=1
#                                               ,dataset = "cgp2016")
# predictedPtype_Cisplatin <- data.frame(predictedPtype_Cisplatin)
# 
# tcga_durg_ic50_res <- predictedPtype_Cisplatin

drugs <- c("Cisplatin","Erlotinib","Rapamycin","Sunitinib","PHA-665752","MG-132","Paclitaxel","Cyclopamine","AZ628","Sorafenib","VX-680","Imatinib","TAE684","Crizotinib","Saracatinib","S-Trityl-L-cysteine","Z-LLNle-CHO","Dasatinib","GNF-2","CGP-60474","CGP-082996","A-770041","WH-4-023","WZ-1-84","BI-2536","BMS-509744","CMK","Pyrimethamine","JW-7-52-1","A-443654","GW843682X","MS-275","Parthenolide","KIN001-135","TGX221","Bortezomib","XMD8-85","Roscovitine","Salubrinal","Lapatinib","Vinorelbine","NSC-87877","QS11","CP466722","Midostaurin","Shikonin","AKT inhibitor VIII","Embelin","Bexarotene","Bleomycin","Phenformin")
length(drugs)
# for (drug in drugs) {
#   print(drug)
#   set.seed(12345)
#   tmpic50 <- pRRopheticPredict(as.matrix(tcga_exp)
#                                , drug
#                                , selection=1
#                                , dataset = "cgp2016")
#   tmpic50 <- data.frame(tmpic50)
#   colnames(tmpic50) <- drug
#   tcga_durg_ic50_res <- cbind(tcga_durg_ic50_res, tmpic50)
# }
 save(tcga_durg_ic50_res,file='07_drug/tcga_durg_ic50.RData')
load("/home/pub252/users/liy/20250220_UVM/07_drug/tcga_durg_ic50.RData")

tcga_durg_ic50_res[1:5,1:5]
#tcga_durg_ic50_res <- tcga_durg_ic50_res[,-1]
#colnames(tcga_durg_ic50_res)[1]='Cisplatin'

diff_drug<-function(dat,group){
  dat=data.frame(cluster=group,(dat))
  #gr=c('High','Low')
  gr=names(table(group))
  dat1=dat[dat$cluster==gr[1],-1]
  dat2=dat[dat$cluster==gr[2],-1]
  #dat3=dat[dat$cluster==gr[3],-1]
  drug=unique(colnames(dat)[-1])
  p_vale=data.frame()
  for (i in drug){
    # x=c(dat1[,i],dat2[,i],dat3[,i])
    # g= factor(rep(names(table(group)), c(nrow(dat1), nrow(dat2), nrow(dat3))),
    #           labels = names(table(group)))
    dd1=wilcox.test(dat1[,i],dat2[,i])$p.value
    #dd1=kruskal.test(x,g)$p.value
    p_vale=rbind(p_vale,data.frame(drug=i,p.value=dd1))
  }
  return(p_vale)
}
COMMSAMPLE <- intersect(tcga.risktype.cli$Samples,rownames(tcga_durg_ic50_res))
tcga.pahtway.diff<-diff_drug(dat=tcga_durg_ic50_res[COMMSAMPLE,],group=tcga.risktype.cli[COMMSAMPLE,]$Risktype)
rownames(tcga.pahtway.diff)=tcga.pahtway.diff$drug
head(tcga.pahtway.diff)
table(tcga.pahtway.diff$p.value<0.05)
tcga.pahtway.diff <- tcga.pahtway.diff[tcga.pahtway.diff$p.value<0.05,]
drug <- gsub("\\.", "-", rownames(tcga.pahtway.diff))
mg_PlotMutiBoxplot(data = tcga_durg_ic50_res[,rownames(tcga.pahtway.diff)],group = tcga.risktype.cli$Risktype,
                   legend.pos = 'top',group_cols = risktype.col,add = 'boxplot',test_method = 'wilcox.test',ylab = 'IC50')
tcga.drug.dat=cbind(tcga.risktype.cli[,'Riskscore'],tcga_durg_ic50_res[tcga.risktype.cli$Samples,])
drug [14]<- "Obatoclax Mesylate"
tcga_tide_list <- list()
for (fea in drug) {
  #fea <- drug[2]
  print(fea)
  tmp_plot <- mg_violin_1(data.frame(tcga.risktype.cli[COMMSAMPLE,]$Risktype
                                     ,tcga_durg_ic50_res[COMMSAMPLE, fea])
                          ,melt = T
                          ,ylab = fea
                          ,jitter=T
                          ,group_col = risktype.col#pal_jco()(9)[1:2]
                          ,test_method = 'wilcox.test'
                          ,cmp_test_method = 'wilcox.test'
                          
                          ,show_compare = T,
                          legend.pos=NULL)
  tcga_tide_list[[fea]] <- tmp_plot
}
fig8<- cowplot::plot_grid(plotlist = tcga_tide_list,
                            ncol = 6,nrow = 3)
fig8

ggsave("07_drug/fig7.pdf",fig8,height = 15,width =18)





