---
title: "Dolen_OT_neurons_enrichment_tests"
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
library("GSEABase")
#library(topGO)
#library(ALL)
library("readxl")
#data(ALL)
#data(geneList)
```


```{r objects}
dat<-readRDS("../data/filtered10kOxy.rds")
pData(dat)$CellType<-"Magnocellular"
pData(dat)$CellType[pData(dat)$Cluster==1]<-"Parvocellular"
PV_vs_Magno_res<-readRDS("cache/PV_vs_Magno_res.rds")
nCellCutoff<-5
expressed_genes<-rownames(fData(dat))[fData(dat)$num_cells_expressed>=nCellCutoff]

qval_cutoff<-0.001
PV_vs_Magno_sigGeneIDs<-PV_vs_Magno_res$gene_id[PV_vs_Magno_res$qval<=qval_cutoff]
length(PV_vs_Magno_sigGeneIDs)

#Calculate Cell Type Specific mean expression
fData(dat)$Magno_mean_cpc<-Matrix::rowMeans(exprs(dat)[,pData(dat)$CellType=="Magnocellular"])
fData(dat)$Parvo_mean_cpc<-Matrix::rowMeans(exprs(dat)[,pData(dat)$CellType=="Parvocellular"])
fData(dat)$logfc<-log2(fData(dat)$Magno_mean_cpc/fData(dat)$Parvo_mean_cpc)


#Cluster_res<-differentialGeneTest(dat[expressed_genes,],fullModelFormulaStr = "~num_genes_expressed+Cluster",reducedModelFormulaStr = "~num_genes_expressed",cores=6)
#qval_cutoff<-0.001
#Cluster_sigGeneIDs<-Cluster_res$gene_id[Cluster_res$qval<=qval_cutoff]
#length(Cluster_sigGeneIDs)
#magExp<-apply(exprs(dat)[Cluster_sigGeneIDs,pData(dat)$CellType=="Magnocellular"],1,mean)
#pavExp<-apply(exprs(dat)[Cluster_sigGeneIDs,pData(dat)$CellType=="Parvocellular"],1,mean)
#delta<-magExp-pavExp
#CT<-ifelse(delta<0,"Parvo","Magno")


```


```{r gene_lists}



###################################################################################
#genes associated with schizophrenia with p<5-08 and p<1e-06 in a meta-analysis of PGC GWAS datasets
schiz <- read_excel("lists4enrichment/Schizophrenia_2018_21_MOESM4_ESM.xlsx", sheet = "a.PRS1andPRS2genes(distance)")
schiz$SYMBOL

#use homologene
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

map<-getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = schiz$SYMBOL, mart = human,
      attributesL = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position"), martL = mouse)
str(map)

indx<-match(schiz$SYMBOL,map$HGNC.symbol )
schiz$SYMBOL[1]; map$HGNC.symbol[indx[1]]
indx<-indx[!is.na(indx)]
map<-map[indx,]             
map<-map[map$MGI.symbol %in% lookupGeneName(dat,expressed_genes),]
str(map)

#test enrichment 
GeneSetToTest<-map$MGI.symbol
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)

#Are there more GeneSetToTest in parvo relative magno as compared to the distribution of all expressed genes
wilcox.test(fData(dat[lookupGeneId(dat, GeneSetToTest),])$logfc, 
	fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)), ])$logfc,
	 alternative = "less")


#Are there more intersect(GeneSetToTest,significant_DE_genes) in parvo relative magno as 
#compared to the distribution of all expressed genes

#wilcox.test(fData(dat[lookupGeneId(dat,intersect(GeneSetToTest,
#	lookupGeneName(dat, PV_vs_Magno_sigGeneIDs))), ])$logfc, 
#    fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)),])$logfc, 
#    alternative = "less")

schiz<-GeneSetToTest

###################################################################################
#type 2 diabetes, alzheimers

t2d<-read.csv("lists4enrichment/DiabetesType2.txt", stringsAsFactors=FALSE,skip=1)
str(t2d)

map<-getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = t2d[,1], mart = human,
      attributesL = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position"), martL = mouse)
str(map)

indx<-match(t2d[,1],map$HGNC.symbol )
map<-map[indx,]             
map<-map[map$MGI.symbol %in% lookupGeneName(dat,expressed_genes),]

#test enrichment 
GeneSetToTest<-map$MGI.symbol
#GeneSetToTest<-t2d
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)

#Are there more GeneSetToTest in parvo relative magno as compared to the distribution of all expressed genes
wilcox.test(fData(dat[lookupGeneId(dat, GeneSetToTest),])$logfc, 
	fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)), ])$logfc,
	 alternative = "less")

#Are there more intersect(GeneSetToTest,significant_DE_genes) in parvo relative magno as 
#compared to the distribution of all expressed genes
wilcox.test(fData(dat[lookupGeneId(dat,intersect(GeneSetToTest,
	lookupGeneName(dat, PV_vs_Magno_sigGeneIDs))), ])$logfc, 
    fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)),])$logfc, 
    alternative = "less")


t2d<-GeneSetToTest

###################################################################################
# alzheimers

alz<-read.delim("lists4enrichment/Alzheimers-gwas-association-downloaded_2019-10-21-EFO_0006514-withChildTraits.tsv",stringsAsFactors=FALSE)
str(alz)
alz_genes<-unique(alz$REPORTED.GENE.S.)

#remove multiple entries per snp and duplicates
alz_genes<-c(alz_genes,unique(unlist(sapply(alz_genes[grep(",",alz_genes)],strsplit,","))))
alz_genes<-alz_genes[-grep(",",alz_genes)]
grep(",",alz_genes)

#use homologene
human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")

map<-getLDS(attributes = c("ensembl_gene_id","hgnc_symbol"),filters = "hgnc_symbol", values = alz_genes, mart = human,
      attributesL = c("ensembl_gene_id","mgi_symbol","chromosome_name","start_position"), martL = mouse)
str(map)

indx<-match(alz_genes,map$HGNC.symbol )
map<-map[indx,]             
map<-map[map$MGI.symbol %in% lookupGeneName(dat,expressed_genes),]


#test enrichment 
GeneSetToTest<-map$MGI.symbol
#GeneSetToTest<-alz_genes
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)

#Are there more GeneSetToTest in parvo relative magno as compared to the distribution of all expressed genes
wilcox.test(fData(dat[lookupGeneId(dat, GeneSetToTest),])$logfc, 
	fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)), ])$logfc,
	 alternative = "less")

#Are there more intersect(GeneSetToTest,significant_DE_genes) in parvo relative magno as 
#compared to the distribution of all expressed genes

wilcox.test(fData(dat[lookupGeneId(dat,intersect(GeneSetToTest,
	lookupGeneName(dat, PV_vs_Magno_sigGeneIDs))), ])$logfc, 
    fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)),])$logfc, 
    alternative = "less")


alz_genes<-GeneSetToTest


###################################################################################
# Updated SAFARI

usgl<-read.csv("../data/4enrichment/SFARI-Gene_genes_05-19-2020release_06-16-2020export.csv",stringsAsFactors=FALSE)
str(usgl)

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

#test enrichment 
GeneSetToTest<-map$MGI.symbol
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)

#Are there more GeneSetToTest in parvo relative magno as compared to the distribution of all expressed genes
wilcox.test(fData(dat[lookupGeneId(dat, GeneSetToTest),])$logfc, 
  fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)), ])$logfc,
   alternative = "less")

#Are there more intersect(GeneSetToTest,significant_DE_genes) in parvo relative magno as 
#compared to the distribution of all expressed genes

wilcox.test(fData(dat[lookupGeneId(dat,intersect(GeneSetToTest,
  lookupGeneName(dat, PV_vs_Magno_sigGeneIDs))), ])$logfc, 
    fData(dat[lookupGeneId(dat, lookupGeneName(dat, expressed_genes)),])$logfc, 
    alternative = "less")

usgl<-GeneSetToTest


###################################################################################
#ASDRisk
NewASDRiskGenes<-readxl::read_excel("../annotations/ASDRiskGenes_Ale_update_6-5-19.xlsx",sheet="ASD Risk Genes",skip=1)
NewASDRiskGenes$Gene<-tools::toTitleCase(stringr::str_to_lower(NewASDRiskGenes$Gene))
NewASDRiskGenes<-NewASDRiskGenes[NewASDRiskGenes$Gene %in% lookupGeneName(dat,expressed_genes),]
categoryCutoff<-4
NewASDRiskGenes<-NewASDRiskGenes[NewASDRiskGenes$Category<=categoryCutoff,]
NewASDRiskGenes

#ASDRisk
GeneSetToTest<-NewASDRiskGenes$Gene
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))

phyper(q,m,n,k, lower.tail=FALSE)


###################################################################################
#FMRP list 
FMRPgenes<-readxl::read_excel("../annotations/FMRP binding partners Darnell.xlsx")
FMRPgenes<-FMRPgenes[FMRPgenes$'gene symbol' %in% lookupGeneName(dat,expressed_genes),]
FMRPgenes

#FMRP
GeneSetToTest<-FMRPgenes$Gene
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)

# Intersection
GeneSetToTest<-intersect(FMRPgenes$Gene, NewASDRiskGenes$Gene)
k<-length(PV_vs_Magno_sigGeneIDs)
n<-length(expressed_genes)-length(GeneSetToTest)
m<-length(GeneSetToTest)
q<-length(intersect(GeneSetToTest,lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)))
phyper(q,m,n,k, lower.tail=FALSE)





```

```{r plot_densities}

#remove genes with logfc == Inf
fData(dat)<-fData(dat)[-which(fData(dat)$logfc==Inf),]


fData(dat)$alz_genes<-fData(dat)$gene_short_name %in% alz_genes
fData(dat)$t2d<-fData(dat)$gene_short_name %in% t2d
fData(dat)$schiz<-fData(dat)$gene_short_name %in% schiz
fData(dat)$isSignificant<-fData(dat)$gene_short_name %in% lookupGeneName(dat,PV_vs_Magno_sigGeneIDs)


p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") + 
  geom_density(aes(x=logfc),color="green",fill="green",alpha=0.3,data=fData(dat)[fData(dat)$schiz,]) + 
  geom_density(aes(x=logfc),color="darkgreen",fill="darkgreen",alpha=0.3,data=fData(dat)[fData(dat)$schiz & fData(dat)$isSignificant,]) +
  geom_vline(aes(xintercept=0),linetype="dashed") +
  xlim(-10,10)+ ggtitle("Density of Schizophrenia Risk Genes") +
  geom_rug(aes(x=logfc),color="darkgreen",fill="darkgreen",alpha=0.3,data=fData(dat)[fData(dat)$schiz & fData(dat)$isSignificant,])
p

pdf("../plots/logfc_distribution_by_SchizophreniaRiskGene.pdf",width=5,height=5)
p
dev.off()

p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") +
  geom_density(aes(x=logfc),color="hotpink",fill="hotpink",alpha=0.3,data=fData(dat)[fData(dat)$t2d,]) + 
  geom_density(aes(x=logfc),color="maroon2",fill="maroon2",alpha=0.3,data=fData(dat)[fData(dat)$t2d & fData(dat)$isSignificant,]) + 
   xlim(-10,10)+  ggtitle("Density of Type II Diabetes Genes") +
  geom_vline(aes(xintercept=0),linetype="dashed") +
  geom_rug(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$t2d & fData(dat)$isSignificant,])
  monocle:::monocle_theme_opts()
p

pdf("../plots/logfc_distribution_by_TypeIIDiabetes.pdf",width=5,height=5)
p
dev.off()

fData(dat)[fData(dat)$t2d[rank(fData(dat)[fData(dat)$t2d,"logfc"])],]


p<-ggplot(fData(dat[expressed_genes,])) +
  geom_density(aes(x=logfc),alpha=0,color="black") + 
  geom_density(aes(x=logfc),color="gold",fill="gold",alpha=0.3,data=fData(dat)[fData(dat)$alz_genes,]) + 
  geom_density(aes(x=logfc),color="goldenrod2",fill="goldenrod2",alpha=0.3,data=fData(dat)[fData(dat)$alz_genes & fData(dat)$isSignificant,]) +
  geom_vline(aes(xintercept=0),linetype="dashed") +
  xlim(-10,10)+ ggtitle("Density of Alzheimers Risk Genes") +
  geom_rug(aes(x=logfc),color="cyan",fill="cyan",alpha=0.3,data=fData(dat)[fData(dat)$alz_genes & fData(dat)$isSignificant,])
  monocle:::monocle_theme_opts()
p

pdf("../plots/logfc_distribution_by_AlzheimersRiskGene.pdf",width=5,height=5)
p
dev.off()



```




# Session Information
```{r}
sessionInfo()
```