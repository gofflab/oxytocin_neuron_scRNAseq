#figs for manuscript
---
title: "Dolen_OT_neurons_figs"
author: "Genevieve Stein-O'Brien"
date: "06/18/2020"
output: html_document
---

# Initialize Environment
```{r init}
source('code/init.R')
library(reticulate)
use_condaenv('OT-reticulate')
set.seed(1234)
library("AnnotationDbi")
library("org.Mm.eg.db")
library("biomaRt")
library("RColorBrewer")
library(corrplot)
library(monocle)
library("yarrr")
```

```{r load_data}
dat<-readRDS(file="../data/filtered10kOxy.rds")

#call cell types by clusters
pData(dat)$CellType<-"Magnocellular"
pData(dat)$CellType[pData(dat)$Cluster==1]<-"Parvocellular"
PMcolors<-c("orangered","darkcyan")

PV_vs_Magno_res<-readRDS("PV_vs_Magno_res.rds")
nCellCutoff<-5
expressed_genes<-rownames(fData(dat))[fData(dat)$num_cells_expressed>=nCellCutoff]

#qval_cutoff<-0.0000001
qval_cutoff<-0.001
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
length(PV_vs_Magno_sigGeneIDs)
```

```{r fig5bc}

grid.arrange(
	plot_cell_clusters(dat,color="CellType", shape="source_plate") + scale_color_manual(values=PMcolors),
	plot_cell_clusters(dat,color="Fluorogold", shape="source_plate") + scale_color_manual(values=label_colors),
	ncol=2,nrow=1, widths = c(1, 1), heights=c(1))
)


   p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos)) + 
   	  geom_point(aes_string(color=color_by,alpha="value"),stroke=0,size=cell_size) + 
      facet_wrap('gene_short_name')+ theme_bw() + scale_color_manual(values=PMcolors) + 
      monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 



pData(dat)$Fluorogold
pData(dat)$Cluster
sort_date
source_plate


 
#coordinates for plotting 
dat@reducedDimS


```



```{r variance_test}

parvo.dis<-dist(t(exprs(dat[expressed_genes,pData(dat)$CellType=='Parvocellular'])))
parvo.distData<-data.frame(celltype="Parvocellular",JSD=as.vector(parvo.dis))

magno.dis<-dist(t(exprs(dat[expressed_genes,pData(dat)$CellType=='Magnocellular'])))
magno.distData<-data.frame(celltype="Magnocellular",JSD=as.vector(magno.dis))

#rand.dis<-dist(t(exprs(dat[expressed_genes,rownames(pData(dat)) %in% sample(rownames(pData(dat)),30)])))

rand.dis<-vector()
for(i in 1:10){
temp<-dist(t(exprs(dat[expressed_genes,rownames(pData(dat)) %in% sample(rownames(pData(dat)),21)])))
rand.dis<-c(rand.dis,temp)
}
#rand.dis<-dist(t(exprs(dat[expressed_genes,rownames(pData(dat))])))
rand.distData<-data.frame(celltype="rand",JSD=as.vector(rand.dis))
rand.distData.ecdf<-ecdf(rand.distData$JSD)

distData<-rbind(parvo.distData,magno.distData,rand.distData)


# DKW confidence intervals
n<-200
alpha<-0.05
eps<-sqrt(log(2/alpha)/(2*n))
#xx            <-    seq(60,120,length.out=1000)
xx      <-  knots(rand.distData.ecdf)
ll            <-    pmax(rand.distData.ecdf(xx)-eps,0)
uu             <-    pmin(rand.distData.ecdf(xx)+eps,1)
randCI<-as.data.frame(cbind(xx,ll,uu))

# ks test
ks.test(parvo.distData$JSD,magno.distData$JSD,exact=T)

# Difference in means
t.test(parvo.distData$JSD,magno.distData$JSD)

#Variance test
var.test(parvo.distData$JSD,magno.distData$JSD,alternative="greater")

distData<-rbind(parvo.distData,magno.distData)

p<-ggplot(distData) + 
  #geom_ribbon(aes(x=xx,ymin=ll,ymax=uu),color="grey80",fill="grey80",alpha=0.25,data=randCI) + 
  stat_ecdf(aes(x=JSD,color=celltype),geom="line") + 
  monocle:::monocle_theme_opts() + 
  scale_color_manual(values=c(PMcolors,"black")) + 
  ggtitle("Intra-celltype heterogeneity of gene expression") + 
  ylab("ecdf") + xlab("Euclidean Distance") +
  geom_hline(yintercept=1,linetype="dashed",color="grey70") + 
  geom_hline(yintercept=0,linetype="dashed",color="grey70") +
  theme(legend.position = c(0.8, 0.2))
p

pdf("../plots/ecdf.pdf",width=5,height=5)
p + theme(legend.position = c(0.8, 0.5))
dev.off()

ks.test(distData$JSD[distData$celltype==0],distData$JSD[distData$celltype==1])

```





```{r fig5bc}

genes2plot <- dat[row.names(subset(fData(dat), gene_short_name %in% c("Kcnmb4","Calb1","Cnr1","Reln"))),]

pdf("../plots/fig5bc.pdf",width=3.17,height=16,fonts="Helvetica")
grid.arrange(
	plot_cell_clusters(dat,color="CellType",size=3) + scale_color_manual(values=PMcolors),
	plot_cell_clusters(dat,color="Fluorogold") + scale_color_manual(values=label_colors),
	plot_cell_clusters(dat,markers=c("Kcnmb4","Calb1","Cnr1","Reln")) + 
		scale_color_viridis(option="viridis") +scale_size_area(max_size = 6),
	monocle::plot_genes_violin(genes2plot, grouping="CellType", ncol=2, min_expr=0.1,
		color_by="CellType",panel_order=c("Kcnmb4","Calb1","Cnr1","Reln")) +
		scale_fill_manual(values=PMcolors) + theme(legend.position = "top"),
	ncol=2,nrow=3, widths = c(1, 1), heights=c(1,2,2),
  	layout_matrix = rbind(c(1,2),c(3,3),c(4,4))
)
dev.off()
```

```{r fig3a}

myPHeatmap<-function(cds,geneset,logMode=TRUE){
  sub<-cds[lookupGeneId(cds,geneset),]
  mat<-as.matrix(exprs(sub))
  if(logMode){
    mat<-log10(mat+1)
  }
  pheatmap(mat=log10(as.matrix(exprs(sub))+1),
            scale="row",
            labels_row=fData(sub)$gene_short_name,
            annotation_col=pData(sub)[,c("source_plate","Total_mRNAs","num_genes_expressed","Cluster","Fluorogold")],
            labelsCol=FALSE,  
            #annotation_colors=c(label_colors,PMcolors),
            color = colorRampPalette(piratepal(palette="brave"))(100), #rev(brewer.pal(n = 11, name ="RdBu")))(100),
            #color = magma(100),
            clustering_distance_cols = "correlation",
            clustering_distance_rows = "correlation",
            show_rownames = F, show_colnames = F,
            cutree_rows=2, cutree_cols=2, breaks=seq(-3,3,length=101))
}

HM<-myPHeatmap(dat,PV_vs_Magno_sigGeneIDs,logMode=TRUE)
```

```{r density_plots}

ASDRiskGenes<-readxl::read_excel("../annotations/ASDRiskGenes.xlsx",sheet="FMRP ASD Risk Genes",skip=1)
ASDRiskGenes
ASDRiskGenes<-ASDRiskGenes[ASDRiskGenes$Gene %in% lookupGeneName(dat,expressed_genes),]

SFARIRiskGenes<-readxl::read_excel("../annotations/ASDRiskGenes.xlsx",sheet="SFARI ASD RISK GENES",skip=2)
SFARIRiskGenes
SFARIRiskGenes$mouseGeneName<-tools::toTitleCase(stringr::str_to_lower(SFARIRiskGenes$Gene))
SFARIRiskGenes<-SFARIRiskGenes[SFARIRiskGenes$mouseGeneName %in% lookupGeneName(dat,expressed_genes),]

FMRPgenes<-readxl::read_excel("../annotations/FMRP binding partners Darnell.xlsx")
FMRPgenes
FMRPgenes<-FMRPgenes[FMRPgenes$'gene symbol' %in% lookupGeneName(dat,expressed_genes),]

#use homologene for SFARIR
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

map<-getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = SFARIRiskGenes$Gene, mart = human,
      attributesL = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position"), martL = mouse)

#Calculate Cell Type Specific mean expression
fData(dat)$Magno_mean_cpc<-Matrix::rowMeans(exprs(dat)[,pData(dat)$CellType=="Magnocellular"])
fData(dat)$Parvo_mean_cpc<-Matrix::rowMeans(exprs(dat)[,pData(dat)$CellType=="Parvocellular"])
fData(dat)$logfc<-log2(fData(dat)$Magno_mean_cpc/fData(dat)$Parvo_mean_cpc)

fData(dat)$logfc[c("Calb1", "Kcnmb4", "Reln", "Cnr1")]

fData(dat)$Magno_mean_cpc[c("Calb1", "Kcnmb4", "Reln", "Cnr1")]

qval_cutoff<-0.1
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]


fData(dat[lookupGeneId(dat,intersect(ASDRiskGenes$Gene,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))),])$logfc
fData(dat[lookupGeneId(dat,intersect(map$MGI.symbol,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))),])$logfc
fData(dat[lookupGeneId(dat,intersect(FMRPgenes$'gene symbol',lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))),])$logfc

fData(dat)$ASDRiskGene<-fData(dat)$gene_short_name %in% intersect(FMRPgenes$'gene symbol', map$MGI.symbol)
fData(dat)$SFARIRiskGene<-fData(dat)$gene_short_name %in% SFARIRiskGenes$mouseGeneName
fData(dat)$SFARIRiskGene<-fData(dat)$gene_short_name %in% map$MGI.symbol
fData(dat)$FMRPgenes<-fData(dat)$gene_short_name %in% FMRPgenes$'gene symbol'
fData(dat)$isSignificant<-fData(dat)$gene_short_name %in% lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)



p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") + 
  geom_density(aes(x=logfc),color="orange",fill="orange",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene,]) + 
  geom_density(aes(x=logfc),color="red",fill="red",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant,]) + 
  geom_density(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes,]) + 
  geom_density(aes(x=logfc),color="blue",fill="blue",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes & fData(dat)$isSignificant,]) + 
  geom_density(aes(x=logfc),color="green",fill="green",alpha=0.3,data=fData(dat)[fData(dat)$SFARIRiskGene,]) + 
  geom_density(aes(x=logfc),color="purple",fill="purple",alpha=0.3,data=fData(dat)[fData(dat)$SFARIRiskGene & fData(dat)$isSignificant,]) + 
  geom_point(aes(x=logfc,y=0),color="red",data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant,]) + 
  geom_label_repel(aes(x=logfc,y=0,label=gene_short_name),color="black",data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant,]) + 
  geom_rug(aes(x=logfc),color="orange",fill="orange",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene,]) + 
  geom_rug(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes,]) +
  xlim(-10,10)+
  geom_vline(aes(xintercept=0),linetype="dashed") +
  monocle:::monocle_theme_opts()
p

pdf("../plots/logfc_distribution_by_genessets.pdf",width=5,height=5)
p
dev.off()

Data(dat)[fData(dat)$FMRPgenes & fData(dat)$isSignificant,],binwidth=.5)

p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") + 
  geom_density(aes(x=logfc),color="orange",fill="orange",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene,]) + 
  geom_density(aes(x=logfc),color="red",fill="red",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant,]) +
  geom_vline(aes(xintercept=0),linetype="dashed") +
#  geom_rug(aes(x=logfc),color="orange",fill="orange",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene,]) + 
  geom_rug(aes(x=logfc),color="red",fill="red",alpha=0.3,data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant,]) +
  geom_label_repel(aes(x=logfc,y=0,label=gene_short_name),color="black",data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$isSignificant & fData(dat)$FMRPgenes,]) + 
  xlim(-10,10)+ ggtitle("Density of ASD Risk Genes") +
  monocle:::monocle_theme_opts() 
p

pdf(paste0("../plots/logfc_distribution_by_ASDRiskGene_FDR",qval_cutoff,".pdf"),width=5,height=3)
p
dev.off()


p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") +
  geom_density(aes(x=logfc,),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes,]) + 
  geom_density(aes(x=logfc),color="blue",fill="blue",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes & fData(dat)$isSignificant,],binwidth=.5) + 
#  geom_rug(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes,]) + 
  geom_rug(aes(x=logfc),color="blue",fill="blue",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes & fData(dat)$isSignificant,]) +
  xlim(-10,10)+  ggtitle("Density of FMRP Genes") +
  geom_vline(aes(xintercept=0),linetype="dashed") +
  geom_label_repel(aes(x=logfc,y=0,label=gene_short_name),color="black",data=fData(dat)[fData(dat)$ASDRiskGene & fData(dat)$FMRPgenes & fData(dat)$isSignificant,]) + 
  monocle:::monocle_theme_opts()
p

pdf(paste0("../plots/logfc_distribution_by_FMRP",qval_cutoff,".pdf"),width=5,height=3)
p
dev.off()

p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") +
  geom_density(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes,]) + 
  geom_density(aes(x=logfc),color="blue",fill="blue",alpha=0.3,data=fData(dat)[fData(dat)$FMRPgenes & fData(dat)$isSignificant,]) + 
  #geom_jitter(aes(x=logfc,y=0.3),height=0.02,color="green",data=fData(dat)[fData(dat)$SFARIRiskGene,]) + 
  geom_density(aes(x=logfc),color="green",fill="green",alpha=0.3,data=fData(dat)[fData(dat)$SFARIRiskGene,]) + 
  #geom_jitter(aes(x=logfc,y=0.1),height=0.02,color="red",data=fData(dat)[fData(dat)$SFARIRiskGene & fData(dat)$isSignificant,]) + 
  geom_density(aes(x=logfc),color="purple",fill="purple",alpha=0.3,data=fData(dat)[fData(dat)$SFARIRiskGene & fData(dat)$isSignificant,]) + 
  xlim(-10,10)+
  geom_vline(aes(xintercept=0),linetype="dashed") +
  monocle:::monocle_theme_opts()
p

pdf("../plots/logfc_distribution_by_genessets.pdf",width=5,height=5)
p
dev.off()

```


```{r final_enrichment_lists4sup_tables}

dat<-readRDS("../data/filtered10kOxy.rds")
pData(dat)$CellType<-"Magnocellular"
pData(dat)$CellType[pData(dat)$Cluster==1]<-"Parvocellular"
PV_vs_Magno_res<-readRDS("PV_vs_Magno_res.rds")
nCellCutoff<-5
expressed_genes<-rownames(fData(dat))[fData(dat)$num_cells_expressed>=nCellCutoff]
Cluster_res<-differentialGeneTest(dat[expressed_genes,],fullModelFormulaStr = "~num_genes_expressed+Cluster",reducedModelFormulaStr = "~num_genes_expressed",cores=6)
saveRDS(Cluster_res,"cache/Cluster_res.rds")
Cluster_res<-readRDS("cache/Cluster_res.rds")


### SUP1 differentially expressed genes
qval_cutoff<-0.001
Cluster_sigGeneIDs<-Cluster_res$gene_id[Cluster_res$qval<=qval_cutoff]
length(Cluster_sigGeneIDs)

magExp<-apply(exprs(dat)[Cluster_sigGeneIDs,pData(dat)$CellType=="Magnocellular"],1,mean)
pavExp<-apply(exprs(dat)[Cluster_sigGeneIDs,pData(dat)$CellType=="Parvocellular"],1,mean)

delta<-magExp-pavExp
CT<-ifelse(delta<0,"Parvo","Magno")

SupTable1<-cbind(Cluster_res[Cluster_res$qval<=qval_cutoff,], magExp, pavExp, CT)
write.csv(SupTable1, file=paste0("csvs/Significantly_DE_genes_qval",qval_cutoff,".csv"))


### SUP2 sig ASDrisk genes

#NewASDRiskGenes<-readxl::read_excel("/Users/genesofeve/Desktop/oxytocin_neuron_scRNAseq/annotations/ASDRiskGenes_Ale_update_6-5-19.xlsx",sheet="ASD Risk Genes",skip=1)
#NewASDRiskGenes$Gene<-tools::toTitleCase(stringr::str_to_lower(NewASDRiskGenes$Gene))
#NewASDRiskGenes<-NewASDRiskGenes[NewASDRiskGenes$Gene %in% lookupGeneName(dat,expressed_genes),]
NewASDRiskGenes<-read.csv("../data/4enrichment/SFARI-Gene_genes_05-19-2020release_06-16-2020export.csv",stringsAsFactors=FALSE)

#use homologene
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

map<-getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = usgl$gene.symbol, mart = human,
      attributesL = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position"), martL = mouse)
str(map)

indx<-match(usgl$gene.symbol,map$HGNC.symbol )
map<-map[indx,]             
map<-map[map$MGI.symbol %in% lookupGeneName(dat,expressed_genes),]
str(map)

GeneSetToTest<-map$MGI.symbol
#GeneSetToTest<-NewASDRiskGenes$Gene
qval_cutoff<-.1
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
intersect(NewASDRiskGenes$Gene,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))
indx<-intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))

magExp<-apply(exprs(dat)[match(indx,fData(dat)$gene_short_name),pData(dat)$CellType=="Magnocellular"],1,mean)
pavExp<-apply(exprs(dat)[match(indx,fData(dat)$gene_short_name),pData(dat)$CellType=="Parvocellular"],1,mean)
delta<-magExp-pavExp
CT<-ifelse(delta<0,"Parvo","Magno")
names(CT)<-indx
SupTable2<-cbind(Cluster_res[match(indx,Cluster_res$gene_short_name),], magExp, pavExp,CT)

write.csv(SupTable2, file=paste0("csvs/Significantly_DE_ASDgenes_qval",qval_cutoff,".csv"))


### SUP3 sig FMRP genes

FMRPgenes<-readxl::read_excel("/Users/genesofeve/Desktop/oxytocin_neuron_scRNAseq/annotations/ASDRiskGenes_Updated.xlsx",sheet="FMRPBindingPartners")
#FMRPgenes<-readxl::read_excel("../annotations/FMRP binding partners Darnell.xlsx")
FMRPgenes$Gene 
GeneSetToTest<-FMRPgenes$Gene[FMRPgenes$Gene %in% lookupGeneName(dat,expressed_genes)]
qval_cutoff<-.1
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
indx<-intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs))

magExp<-apply(exprs(dat)[match(indx,fData(dat)$gene_short_name),pData(dat)$CellType=="Magnocellular"],1,mean)
pavExp<-apply(exprs(dat)[match(indx,fData(dat)$gene_short_name),pData(dat)$CellType=="Parvocellular"],1,mean)
delta<-magExp-pavExp
CT<-ifelse(delta<0,"Parvo","Magno")
names(CT)<-indx
SupTable3<-cbind(Cluster_res[match(indx,Cluster_res$gene_short_name),], magExp, pavExp,CT)

write.csv(SupTable3, file=paste0("csvs/Significantly_DE_FMRPgenes_qval",qval_cutoff,".csv"))

```


# mitochondrial analysis 
```{r}

mito_genes<-fData(dat)$gene_id[grepl("^mt-",fData(dat)$gene_short_name)]
pData(dat)$mt_reads <- Matrix::colSums(exprs(dat)[mito_genes,])
pData(dat)$total_reads  <- Matrix::colSums(exprs(dat))
pData(dat)$mito_ratio <- pData(dat)$mt_reads/pData(dat)$total_reads
mito_summary <- pData(dat)%>%
  group_by(CellType) %>%
  summarise(mean=mean(mito_ratio),sd=sd(mito_ratio),se=sd/sqrt(n()))

pdf("../plots/mito_ratio_barplot.pdf",width=5,height=10)
p<-ggplot(mito_summary) + 
  geom_bar(aes(y=mean,fill=CellType,x=CellType),stat="identity",position="dodge") + 
  geom_errorbar(aes(y=mean,x=CellType,color=CellType,ymin=mean-se,ymax=mean+se),color="black",width=.75,position="dodge") + 
  #geom_boxplot(aes(x=CellType,y=mito_ratio,fill=CellType),notch=FALSE) +
  geom_point(aes(x=CellType,y=mito_ratio,color=Fluorogold),data=pData(dat),position = "jitter") +
  scale_fill_manual(values=c("steelblue","red")) +
  scale_y_continuous(breaks = seq(0, 0.3, by = 0.025)) + 
  scale_color_manual(values=label_colors) +
  #ylim(c(0,0.15)) +
  ylab("Mitochondrial Ratio") + 
  monocle:::monocle_theme_opts()
p
dev.off()

mito_summary <- pData(dat)%>%
  group_by(CellType,source_plate) %>%
  summarise(mean=mean(mito_ratio),sd=sd(mito_ratio),se=sd/sqrt(n()))

pdf("../plots/mito_ratio_boxplots_by_plate.pdf",width=5,height=5)
p<-ggplot(pData(dat)) + 
  #geom_bar(aes(y=mean,fill=CellType,x=source_plate),position="dodge",stat="identity") + 
  #geom_boxplot(aes(y=mean,x=source_plate,fill=CellType,color=CellType,ymin=mean-se,ymax=mean+se),position="dodge",width=.75) + 
  geom_boxplot(aes(x=source_plate,y=mito_ratio,fill=CellType),notch=FALSE) +
  scale_fill_manual(values=c("steelblue","red")) +
  scale_color_manual(values=c("steelblue","red")) +
  ylab("Mitochondrial Ratio") + 
  monocle:::monocle_theme_opts()
p
dev.off()

t.test(pData(dat)$mito_ratio[pData(dat)$CellType=="Magnocellular"],pData(dat)$mito_ratio[pData(dat)$CellType=="Parvocellular"],alternative="greater",paired = FALSE)

```




```{r ref}
myTSNEPlotAlpha<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=2){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,alpha="value"),stroke=0,size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }else{
      p + geom_point(aes_string(color=color_by,alpha="value",stroke=0,shape=shape_by),size=cell_size) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() + scale_alpha(range=c(0.05,1)) 
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by),size=cell_size) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}
```






