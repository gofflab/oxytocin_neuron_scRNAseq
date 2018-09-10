library(monocle)
library(reshape2)
library(ggplot2)
library(marray)
library(stringr)
library(heatmap.plus)
library(tidyr)
library(pheatmap)
library(Rtsne)
require(Heatplus)
require(Hmisc)
library(RColorBrewer)
library(gplots)
library(MASS)
library(vegan)
library(slackr)
library(ggrepel)
library(ggbiplot)
library(gridExtra)
library(tsne)
library(corrplot)
library(dplyr)
library(mclust)
library(ape)
library(flashClust)
library(WGCNA)
library(viridis)
require(clusterProfiler)
require(org.Mm.eg.db)
require(ReactomePA)
#############
#slackr
#############
library(slackr)
slackrSetup(
 channel="#brown_layer_6",
 incoming_webhook_url="https://hooks.slack.com/services/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",
 username="slackR",
 api_token="xoxp-9317854839-9317673428-12053677172-a0d993cf53",
 icon_emoji = ":computer:",
 config_file="")

############
# Color Palette
############
celltype_colors<-c("chocolate1","blue4")
label_colors<-c("red","darkgreen")

############
# Helper functions
############
save.xlsx <- function (file, x, ...)
{
  require(xlsx, quietly = TRUE)
  
  for (i in 1:length(x)) {
    if (i == 1)
      write.xlsx(x[[i]], file, sheetName = names(x)[i], ...)
    else write.xlsx(x[[i]], file, sheetName = names(x)[i],
                    append = TRUE, ...)
  }
}

save.xlsx.2 <- function (file, x, ...)
{
  require(xlsx, quietly = TRUE)
  
  for (i in 1:length(x)) {
    if (i == 1)
      write.xlsx(x[[i]], file, sheetName = paste("Cluster",i,sep=" "), ...)
    else write.xlsx(x[[i]], file, sheetName = paste("Cluster",i,sep=" "),
                    append = TRUE, ...)
  }
}

hclust2 <- function(x, method="ward", ...)
  hclust(x, method="ward.D", ...)

dist2s <- function(x, ...)
  as.dist(1-cor(t(x), method="spearman"))

dist2p <- function(x, ...)
  as.dist(1-cor(t(x), method="pearson"))

dist4e <- function(x, ...)
  dist(1-x, method='euclidean');

dist4m <- function(x, ...)
  dist(1-x, method='manhattan');


myHeatCols<-maPalette(low="steelblue",mid="white",high="darkred",k=100)

myCoolerCols<-maPalette(low="black",mid="darkgreen",high="white",k=100)
#myCoolerCols<-maPalette(low="white",mid="red",high="darkred",k=100)

standardize <- function(z) {
  rowmed <- apply(z, 1, median)
  rowmad <- apply(z, 1, mad)  # median absolute deviation
  rv <- sweep(z, 1, rowmed,"-")  #subtracting median expression
  rv <- sweep(rv, 1, rowmad, "/")  # dividing by median absolute deviation
  return(rv)
}

lookupGeneId<-function(eset,gene_names){
  res <- rownames(fData(eset))[fData(eset)$gene_short_name %in% gene_names]
  res <- c(res,rownames(fData(eset))[rownames(fData(eset)) %in% gene_names])
  res <- unique(res)
  res
}

lookupGeneName<-function(eset,gene_id){
  res <- fData(eset[gene_id,])$gene_short_name
  res <- unique(res)
  res
}

#source('hNSC_funcs.R')

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

# myHeatmap<-function(cds,geneset,pseudotime=FALSE){
#   sub<-cds[lookupGeneId(cds,geneset),]
#   if(pseudotime){
#     heatmap.2(log10(exprs(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])+1),scale="none",trace="none",col=brewer.pal(9,"GnBu"),ColSideColors=gg_color_hue(5)[pData(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])$State],labRow=fData(sub[,order(pData(sub)$Pseudotime,decreasing=FALSE)])$gene_short_name,Colv=FALSE,labCol=FALSE)#,distfun=function(x){JSdist(t(x))})
#   }else {
#     heatmap.2(log10(exprs(sub)+1),scale="none",trace="none",col=brewer.pal(9,"GnBu"),ColSideColors=gg_color_hue(5)[pData(sub)$State],labRow=fData(sub)$gene_short_name,labCol=FALSE)#,distfun=function(x){JSdist(t(x))})
#   }
# }

myHeatmap<-function(cds,geneset,logMode=FALSE){
  sub<-cds[lookupGeneId(cds,geneset),]
  if(logMode){
    heatmap.2(log10(exprs(sub)+1),scale="none",trace="none",col=myCoolerCols,labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(3)[pData(sub)$kmeans_tSNE_cluster],labCol=FALSE,distfun=dist,hclustfun=hclust2)
  }else{
    heatmap.2(exprs(sub),scale="none",trace="none",col=myCoolerCols,labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(3)[pData(sub)$kmeans_tSNE_cluster],labCol=FALSE,distfun=dist,hclustfun=hclust2)
  }
}

meltCDS<-function(cds,geneset,logMode=F){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-as.matrix(exprs(sub))
  if(logMode){
    sub.expr<-log10(sub.expr+1)
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","cell_id","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="cell_id",by.y="cell_id")
  res<-merge(res,fData(sub),by.x="gene_id",by.y="gene_id")
  #print(head(res))
  res
}

vstMeltCDS<-function(cds,geneset,scale=FALSE){
  sub<-cds[lookupGeneId(cds,geneset),]
  sub.expr<-as.matrix(vstExprs(sub))
  if(scale=="row"){
    #sub.expr<-log10(sub.expr/rowMeans(sub.expr))
    sub.expr<-t(scale(t(sub.expr)))
  }
  sub.expr.melt<-melt(sub.expr)
  colnames(sub.expr.melt)<-c("gene_id","cell_id","value")
  res<-merge(sub.expr.melt,pData(sub),by.x="cell_id",by.y="cell_id")
  res<-merge(res,fData(sub),by.x="gene_id",by.y="gene_id")
  #print(head(res))
  res
}

myBarMap<-function(cds,geneset,facet_by="cluster",color_by="factor(cluster)",cluster="both",showSummary=T,scale=FALSE,scale_axes="free",space="free_x",...){
  sub.melt<-meltCDS(cds,geneset,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  sub.melt.summary<-sub.melt %>%
    dplyr::group_by_(.dots=c("gene_short_name",facet_by_melt)) %>%
    dplyr::summarise(mean=mean(value),median=median(value),sd=sd(value),upper_bound=mean+sd,lower_bound=max(mean-sd,0))
  
  if(cluster %in% c("row","both",T)){
    sub.sum.mat<-sub.melt.summary %>%
      recast(as.formula(paste("gene_short_name ~",facet_by)),measure.var="mean",fun.aggregate=mean)
    sub.sum.hclust<-hclust2(dist(sub.sum.mat[,-1]))
    gene.order.idx<-order.dendrogram(as.dendrogram(sub.sum.hclust))
    gene.order<-sub.sum.mat$gene_short_name[gene.order.idx]
    sub.melt$gene_short_name<-factor(sub.melt$gene_short_name, levels=gene.order)
  }
  
  if(cluster %in% c("column","both",T)){
    sub.mat<-sub.melt %>%
      recast(as.formula("gene_short_name ~ cell_id"),measure.var="value",fun.aggregate=mean)
      sub.hclust<-hclust2(dist(t(sub.mat[,-1])))
      cell.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
      cell.order<-colnames(sub.mat[,-1])[cell.order.idx]
      #print(cell.order)
      sub.melt$cell_id<-factor(sub.melt$cell_id,levels=cell.order)
  }
  
  p<-ggplot(sub.melt)
  p<-p + geom_bar(aes_string(x="cell_id",y="value",fill=color_by,color=color_by),stat="identity")
  
  if(showSummary){
    p<-p + geom_hline(aes(yintercept=mean),data=sub.melt.summary,size=1.0)
    p<-p + geom_hline(aes(yintercept=upper_bound),data=sub.melt.summary,linetype="dashed")
    p<-p + geom_hline(aes(yintercept=lower_bound),data=sub.melt.summary,linetype="dashed")
  }
  p<-p +
    facet_grid(as.formula(paste("gene_short_name ~", facet_by)),scale=scale_axes,space=space,labeller=labeller(.default=label_both,gene_short_name=label_value)) +
    theme_bw() + guides(color=FALSE) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.margin = unit(0, "lines"),
          panel.grid = element_blank()
          )
  p
}

myFacetedHeatMap<-function(cds,geneset,facet_by="celltype+Day",cluster="both",scale="row",space="free_x",limits=c(-5,5),...){
  sub.melt<-vstMeltCDS(cds,geneset,scale=scale,...)
  facet_by_melt<-strsplit(facet_by,"\\+")[[1]]
  
  sub.mat<-sub.melt %>%
    recast(as.formula("gene_short_name ~ cell_id"),measure.var="value",fun.aggregate=mean)
  #print(apply(sub.mat[-1],1,mean))
  
  if(cluster %in% c("row","both",T)){
    sub.hclust<-hclust2(dist(sub.mat[,-1]))
    gene.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
    gene.order<-sub.mat$gene_short_name[gene.order.idx]
    sub.melt$gene_short_name<-factor(sub.melt$gene_short_name, levels=gene.order)
  }
  
  if(cluster %in% c("column","both",T)){
    sub.hclust<-hclust2(dist(t(sub.mat[,-1])))
    cell.order.idx<-order.dendrogram(as.dendrogram(sub.hclust))
    cell.order<-colnames(sub.mat[,-1])[cell.order.idx]
    #print(cell.order)
    sub.melt$cell_id<-factor(sub.melt$cell_id,levels=cell.order)
  }
  
  p<-ggplot(sub.melt)
  p<-p + geom_tile(aes_string(x="cell_id",y="gene_short_name",fill="value"))
  
  p<-p +
    facet_grid(as.formula(paste(". ~ ", facet_by)),space=space,scale="free",labeller=labeller(.default=label_both,gene_short_name=label_value)) +
    theme_bw() + guides(color=FALSE) + 
    #scale_fill_viridis(name="Normalized Expression", label=comma) +
    scale_fill_gradient2(low="darkblue",mid="white",high="red",limits=limits,oob=squish) +
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.x=element_blank(),
          strip.text.y = element_text(angle=0,hjust=0),
          strip.background = element_blank(),
          panel.background = element_rect(fill="white"),
          panel.margin = unit(0, "lines"),
          panel.grid = element_blank()
    )
  p
}


plot_genes_line<-function (cds_subset, grouping = "State", min_expr = NULL, cell_size = 0.75, 
          nrow = NULL, ncol = 1, panel_order = NULL, color_by = NULL, label_by_short_name = TRUE, relative_expr = TRUE) 
{
  if (cds_subset@expressionFamily@vfamily %in% c("negbinomial", 
                                                 "negbinomial.size")) {
    integer_expression <- TRUE
  }
  else {
    integer_expression <- FALSE
    relative_expr <- TRUE
  }
  if (integer_expression) {
    cds_exprs <- exprs(cds_subset)
    if (relative_expr) {
      if (is.null(sizeFactors(cds_subset))) {
        stop("Error: to call this function with relative_expr=TRUE, you must call estimateSizeFactors() first")
      }
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/sizeFactors(cds_subset))
    }
    cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  }
  else {
    cds_exprs <- exprs(cds_subset)
    cds_exprs <- reshape2::melt(as.matrix(cds_exprs))
  }
  if (is.null(min_expr)) {
    min_expr <- cds_subset@lowerDetectionLimit
  }
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_exprs$expression[cds_exprs$expression < min_expr] <- min_expr
  cds_pData <- pData(cds_subset)
  cds_fData <- fData(cds_subset)
  cds_exprs <- merge(cds_exprs, cds_fData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_pData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- log10(cds_exprs$expression)
  if (label_by_short_name == TRUE) {
    if (is.null(cds_exprs$gene_short_name) == FALSE) {
      cds_exprs$feature_label <- cds_exprs$gene_short_name
      cds_exprs$feature_label[is.na(cds_exprs$feature_label)] <- cds_exprs$f_id
    }
    else {
      cds_exprs$feature_label <- cds_exprs$f_id
    }
  }
  else {
    cds_exprs$feature_label <- cds_exprs$f_id
  }
  if (is.null(panel_order) == FALSE) {
    cds_exprs$feature_label <- factor(cds_exprs$feature_label, 
                                      levels = panel_order)
  }
  q <- ggplot(aes_string(x = grouping, y = "expression"), data = cds_exprs)
  #if (is.null(color_by) == FALSE) {
  #  q <- q + geom_jitter(aes_string(color = color_by), size = I(cell_size))
  #}
  #else {
  #  q <- q + geom_jitter(size = I(cell_size))
  #}
  
  q <- q + stat_summary(aes_string(color = color_by), fun.data = "mean_cl_boot", 
                        size = 0.35)
  q <- q + stat_summary(aes_string(x = grouping, y = "expression", 
                                   color = color_by, group = color_by), fun.data = "mean_cl_boot", 
                        size = 0.35, geom = "line")

  q <- q + scale_y_log10() + facet_wrap(~feature_label, nrow = nrow, 
                                        ncol = ncol, scales = "free_y")
  if (min_expr < 1) {
    q <- q + expand_limits(y = c(min_expr, 1))
  }
  q <- q + ylab("Expression") + xlab(grouping)
  q <- q + monocle:::monocle_theme_opts()
  q
}


diffAUC <- function(x,y) {
  prediction.use=ROCR::prediction(c(x,y),c(rep(1,length(x)),rep(0,length(y))),0:1)
  perf.use=ROCR::performance(prediction.use,"auc")
  auc.use=round(perf.use@y.values[[1]],3)
  return(auc.use)
}

per_gene_AUC<-function(dat1,dat2,geneset){
  myAUC=unlist(lapply(geneset,function(x)diffAUC(as.numeric(dat1[x,]),as.numeric(dat2[x,]))))
  myAUC[is.na(myAUC)]=0
  avg_diff=unlist(lapply(geneset,function(x)(mean(as.numeric(dat1[x,]))-mean(as.numeric(dat2[x,])))))
  toRet=data.frame(cbind(myAUC,avg_diff),row.names=geneset)
  toRet=toRet[rev(order(toRet$myAUC)),]
  return(toRet)
}

ROC_test<-function(cds,cluster_1,cluster_2,geneset,thresh.use=log(2)) {
  gene_ids<-lookupGeneId(cds,geneset)
  dat<-exprs(cds)
  dat.1<-apply(dat[gene_ids,pData(cds)$cluster %in% cluster_1],1,mean)
  dat.2<-apply(dat[gene_ids,pData(cds)$cluster %in% cluster_2],1,mean)
  total.diff<-abs(dat.1-dat.2)
  genes.diff<-names(which(total.diff>thresh.use))
  genes.use<-genes.diff[genes.diff%in%rownames(dat)]
  res<-per_gene_AUC(dat[,pData(cds)$cluster %in% cluster_1],dat[,pData(cds)$cluster %in% cluster_2],genes.use)
  res<-res[rev(order(abs(res$myAUC-0.5))),]
  res$power=abs(res$myAUC-0.5)*2
  return(res)
}


where.expressed<-function(cds,geneId){
  cpc<-round(exprs(cds[lookupGeneId(cds,geneId),]))
  cpc<- as.vector(cpc > 1)
  cpc<- unlist(lapply(cpc,FUN=function(x){if(x){1}else{0}}))
  cpc
}

jaccard_sim <- function(x) {
  # initialize similarity matrix
  m <- matrix(NA, nrow=ncol(x),ncol=ncol(x),dimnames=list(colnames(x),colnames(x)))
  jaccard <- as.data.frame(m)
  
  for(i in 1:ncol(x)) {
    for(j in i:ncol(x)) {
      jaccard[i,j]= length(which(x[,i] & x[,j])) / length(which(x[,i] | x[,j]))
      jaccard[j,i]=jaccard[i,j]        
    }
  }
}

library(Matrix)
jaccard_distance <- function(m) {
  ## common values:
  A = tcrossprod(m)
  ## indexes for non-zero common values
  im = which(A > 0, arr.ind=TRUE)
  ## counts for each row
  b = rowSums(m)
  
  ## only non-zero values of common
  Aim = A[im]
  
  ## Jacard formula: #common / (#i + #j - #common)
  J = sparseMatrix(
    i = im[,1],
    j = im[,2],
    x = Aim / (b[im[,1]] + b[im[,2]] - Aim),
    dims = dim(A)
  )
  
  return( as.dist(1 - J ))
}

theme_change <- theme(
  plot.background = element_blank(),
  panel.grid.minor = element_blank(),
  panel.grid.major = element_blank(),
  panel.background = element_blank(),
  #panel.border = element_blank(),
  axis.line = element_blank(),
  axis.ticks = element_blank(),
  axis.text.x = element_blank()
  #axis.text.y = element_blank()
  #axis.title.x = element_blank(),
  #axis.title.y = element_blank()
)

myggheatmap<-function(cds,geneset=c("CD44","ACHE","ETV4","ISL1","ISL2","NEUROG2","NEUROD1","CHAT","CRIM1","DIAPH3","BCL11B","FEZF2","FEZF1","LHX1","CRYM","NTNG1","CDH13","CDH22","SOX5","BCL6","NETO1","PCP4","ITM2A","POU3F1"),logMode=T,rowK=5){
  sub<-cds[lookupGeneId(cds,geneset),]
  if(logMode){
    exprs(sub)<-log10(exprs(sub)+1.0)
  }
  
  row.hclust<-hclust(dist(exprs(sub)))
  
  row.cluster<-as.data.frame(cutree(row.hclust,k=rowK))
  colnames(row.cluster)<-c("rowCluster")
  row.order<-order.dendrogram(as.dendrogram(row.hclust))
  
  sub.melt<-melt(exprs(sub))
  
  colnames(sub.melt)<-c("gene_id","sample_id","fpkm")
  
  sub.melt<-merge(sub.melt,pData(sub),by="sample_id")
  
  sub.melt<-merge(sub.melt,fData(sub)[,c("gene_short_name","gene_id")],by="gene_id")
  
  sub.melt<-merge(sub.melt,row.cluster,by.x="gene_id",by.y="row.names")
  
  sub.melt$gene_short_name<-factor(sub.melt$gene_short_name,levels=fData(sub[row.order,])$gene_short_name)
  
  sub.melt$sample_id <- factor(sub.melt$sample_id,levels=pData(sub)$sample_id[order(pData(sub)$Pseudotime,decreasing=F)])
  
  p <- ggplot(sub.melt)
  p + geom_tile(aes(x=sample_id,y=gene_short_name,fill=fpkm)) + facet_grid(rowCluster~State,scales="free",space="free") + theme_change + scale_fill_gradient(low="white",high="darkred") + theme(axis.text.y=element_blank())
}

myCorheatmap<-function(cds,logMode=T,method="color",cor.method="pearson",addrect=NULL,order="hclust",hclust.method="ward",...){
  dat<-as.matrix(vstExprs(cds))
  if(logMode){
    dat<-as.matrix(exprs(cds))
    dat<-log10(dat+1)
  }
  #print(head(t(dat)))
  dat.cor<-cor(t(dat),method=cor.method)
  
  rownames(dat.cor)<-lookupGeneName(cds,rownames(dat.cor))
  colnames(dat.cor)<-lookupGeneName(cds,colnames(dat.cor))
  
  corrplot(dat.cor,method=method,hclust.method=hclust.method,order=order,addrect=addrect,...)
}



PCbiplot <- function(PC, x="PC1", y="PC2", color="black", shape=NULL) {
  # PC being a prcomp object
  data <- data.frame(obsnames=row.names(PC$x), PC$x)
  plot <- ggplot(data, aes_string(x=x, y=y)) + geom_point(alpha=.4, size=3)
  plot <- plot + geom_hline(aes(0), size=.2) + geom_vline(aes(0), size=.2)
  datapc <- data.frame(varnames=rownames(PC$rotation), PC$rotation)
  mult <- min(
    (max(data[,y]) - min(data[,y])/(max(datapc[,y])-min(datapc[,y]))),
    (max(data[,x]) - min(data[,x])/(max(datapc[,x])-min(datapc[,x])))
  )
  datapc <- transform(datapc,
                      v1 = .7 * mult * (get(x)),
                      v2 = .7 * mult * (get(y))
  )
  plot <- plot + coord_equal() + geom_point(data=datapc, aes(x=v1, y=v2, label=varnames), size = 5, vjust=1,color="red")
  plot <- plot + geom_segment(data=datapc, aes(x=0, y=0, xend=v1, yend=v2), arrow=arrow(length=unit(0.2,"cm")), alpha=0.75,color="red") 
  plot <- plot + theme_bw()
  plot
}

myTSNEPlot<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE){
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
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(tmp))
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by,size="value")) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1") + monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,size="value",shape=shape_by)) + facet_wrap('gene_short_name')+ theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }else{
    p<-ggplot(tmp,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(aes_string(color=color_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }else{
      p + geom_point(aes_string(color=color_by,shape=shape_by)) + theme_bw() + scale_color_brewer(palette="Set1")+ monocle:::monocle_theme_opts() 
    }
  }
}

myTSNEPlotAlpha<-function(cds,markers=NULL,logMode=T,color_by="color",shape_by=NULL,scaled=FALSE,cell_size=2){
  tmp<-pData(cds)
  if(!is.null(markers)){
    genes<-as.matrix(exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)]))
    if(logMode){
      genes<-log10(genes+1)
    }
    geneMeans<-rowMax(genes)
    #print(geneMeans)
    if(scaled){
      genes<-genes/geneMeans
    }
    genes<-t(genes)
    #print(genes)
    genes<-melt(genes)
    colnames(genes)<-c("cell_id","gene_id","value")
    genes<-merge(genes,fData(cds),by.x="gene_id",by.y="gene_id",all.x=TRUE,sort=FALSE)
    #print(head(genes))
    tmp<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(tmp))
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

myTSNEPlotRainbow<-function(cds,red="Bdnf",green="Fos",blue="empty",logMode=T,shape_by=NULL,scaled=FALSE,cell_size=2,discrete=FALSE){
  tmp<-pData(cds)
  markers<-c(red,green,blue)
  genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
  if(logMode){
    genes<-log10(genes+1)
  }
  geneMeans<-rowMax(genes)
  if(scaled){
    genes<-genes/geneMeans
  }
  genes<-t(genes)
  colnames(genes)<-lookupGeneName(cds,colnames(genes))
  if(discrete){
    genes<-genes/rowSums(genes)
    genes[is.na(genes)]<-0
  }else{
    genes[is.na(genes)]<-1
  }
  genes<-as.data.frame(genes)
  #Map to Rgb
  if(blue=="empty"){
    genes$empty<-0
  }
  genes$plotColor<-as.vector(unlist(apply(genes[,c(red,green,blue)],1,function(x){rgb(x[1],x[2],x[3])})))
  genes$plotColor[rowSums(genes[,markers])==0]<-"#FFFFFF"
  genes$cell_id<-rownames(genes)
  #print(dim(genes))
  genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
  #print(head(genes))
  titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
  p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
  if(is.null(shape_by)){
    p + geom_point(fill="white",color="black",stroke=0.25,size=cell_size) +
      geom_point(color=genes$plotColor,stroke=0,size=cell_size) +
      monocle:::monocle_theme_opts() + ggtitle(titleString) 
  }else{
    p + geom_point(aes_string(shape=shape_by),fill="white",color="black",stroke=0.25,size=cell_size) +
      geom_point(aes_string(shape=shape_by),stroke=0,color=genes$plotColor,size=cell_size) +
      monocle:::monocle_theme_opts() + ggtitle(titleString) 
  }
}

myTSNEPlotRainbow2<-function(cds,red="Bdnf",green="Fos",blue="empty",logMode=T,shape_by=NULL,scaled=FALSE,cell_size=2,discrete=TRUE){
  tmp<-pData(cds)
  markers<-c(red,green,blue)
  genes<-exprs(cds[rownames(fData(cds)) %in% lookupGeneId(cds,markers)])
  if(logMode){
    genes<-log10(genes+1)
  }
  geneMeans<-rowMax(genes)
  if(scaled){
    genes<-genes/geneMeans
  }
  genes<-t(genes)
  colnames(genes)<-lookupGeneName(cds,colnames(genes))
  if(discrete){
    #genes<-genes/rowSums(genes)
    #genes[is.na(genes)]<-0
    
    genes<-as.data.frame(genes)
    #Map to Rgb
    if(blue=="empty"){
      genes$empty<-0
    }
    genes$cell_id<-rownames(genes)
    #print(dim(genes))
    genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(genes))
    titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
    if(is.null(shape_by)){
      p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
      p + geom_point(fill="white",color = "black",stroke=0.25,size=cell_size,shape=21) +
        geom_point(aes_string(alpha=red),color="red",fill="red",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=green),color="green",fill="green",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=blue),color="blue",fill="blue",stroke=0,size=cell_size) +
        scale_alpha(range=c(0,0.5)) + guides(alpha=FALSE) +
        monocle:::monocle_theme_opts() + ggtitle(titleString) 
    }else{
      p<-ggplot(genes,aes_string(x='tSNE1_pos',y='tSNE2_pos',shape=shape_by))
      p + geom_point(aes_string(shape=shape_by),color="white",stroke=0.25,size=cell_size) +
        geom_point(aes_string(alpha=genes$red),color="red",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=genes$green),color="green",stroke=0,size=cell_size) +
        geom_point(aes_string(alpha=genes$blue),color="blue",stroke=0,size=cell_size) +
        scale_alpha(range=c(0,0.5)) + guides(alpha=FALSE) +
        monocle:::monocle_theme_opts() + ggtitle(titleString) 
    }
  }else{
    genes[is.na(genes)]<-1
 
    genes<-as.data.frame(genes)
    #Map to Rgb
    if(blue=="empty"){
      genes$empty<-0
    }
    genes$plotColor<-as.vector(unlist(apply(genes[,c(red,green,blue)],1,function(x){rgb(x[1],x[2],x[3])})))
    genes$plotColor[rowSums(genes[,markers])==0]<-"#FFFFFF"
    genes$cell_id<-rownames(genes)
    #print(dim(genes))
    genes<-merge(tmp,genes,by.x=0,by.y="cell_id")
    #print(head(genes))
    titleString<-paste("Red=",red,": Green=",green,": Blue=",blue,sep="")
    p<-ggplot(genes,aes(x=tSNE1_pos,y=tSNE2_pos))
    if(is.null(shape_by)){
      p + geom_point(fill="white",color="black",stroke=0.25,size=cell_size) +
        geom_point(color=genes$plotColor,stroke=0,size=cell_size) +
        monocle:::monocle_theme_opts() + ggtitle(titleString) 
    }else{
      p + geom_point(aes_string(shape=shape_by),fill="white",color="black",stroke=0.25,size=cell_size) +
        geom_point(aes_string(shape=shape_by),stroke=0,color=genes$plotColor,size=cell_size) +
        monocle:::monocle_theme_opts() + ggtitle(titleString) 
    }
  }
}
##########
# Geneset import
##########
Arking.Neuron<-read.delim("http://arkinglab.org/upload/GeneLists/NEURO.upper",header=F,stringsAsFactors=F)$V1
Arking.Astro<-read.delim("http://arkinglab.org/upload/GeneLists/ASTRO.upper",header=F,stringsAsFactors=F)$V1
Arking.Synaptic<-read.delim("http://arkinglab.org/upload/GeneLists/SynapticProteins.upper",header=F,stringsAsFactors=F)$V1
CellCycleGenes<-c("ACD","ACTR1A","AHCTF1","AKAP9","ALMS1","ANAPC1","ANAPC10","ANAPC11","ANAPC2","ANAPC4","ANAPC5","ANAPC7","APITD1","ATM","ATR","ATRIP","AURKA","AURKB","AZI1","B9D2","BIRC5","BRCA1","BTRC","BUB1","BUB1B","BUB3","CASC5","CCDC99","CCNA1","CCNA2","CCNB1","CCNB2","CCND1","CCND2","CCND3","CCNE1","CCNE2","CCNH","CDC14A","CDC16","CDC20","CDC23","CDC25A","CDC25B","CDC25C","CDC26","CDC26P1","CDC27","CDC45","CDC6","CDC7","CDCA8","CDK1","CDK2","CDK4","CDK5RAP2","CDK6","CDK7","CDKN1A","CDKN1B","CDKN2A","CDKN2B","CDKN2C","CDKN2D","CDT1","CENPA","CENPC1","CENPH","CENPI","CENPJ","CENPK","CENPL","CENPM","CENPN","CENPO","CENPP","CENPQ","CENPT","CEP135","CEP164","CEP192","CEP250","CEP290","CEP41","CEP57","CEP63","CEP70","CEP72","CEP76","CETN2","CHEK1","CHEK2","CKAP5","CKS1B","CLASP1","CLIP1","CNTRL","CSNK1D","CSNK1E","CUL1","DBF4","DCTN1","DCTN2","DCTN3","DHFR","DHFRP1","DIDO1","DKC1","DNA2","DSN1","DYNC1H1","DYNC1I2","DYNLL1","DYRK1A","E2F1","E2F2","E2F3","E2F4","E2F5","ERCC6L","FBXO5","FEN1","FGFR1OP","FKBP6","GINS1","GINS2","GINS4","GMNN","GORASP1","H2AFX","H2AFZ","HAUS2","HDAC1","HIST1H2AB","HIST1H2AC","HIST1H2AD","HIST1H2AE","HIST1H2AJ","HIST1H2BA","HIST1H2BB","HIST1H2BC","HIST1H2BD","HIST1H2BE","HIST1H2BF","HIST1H2BG","HIST1H2BH","HIST1H2BI","HIST1H2BJ","HIST1H2BK","HIST1H2BL","HIST1H2BM","HIST1H2BN","HIST1H2BO","HIST1H4A","HIST1H4B","HIST1H4C","HIST1H4D","HIST1H4E","HIST1H4F","HIST1H4H","HIST1H4I","HIST1H4J","HIST1H4K","HIST1H4L","HIST2H2AA3","HIST2H2AA4","HIST2H2AC","HIST2H2BE","HIST2H4A","HIST2H4B","HIST3H2BB","HIST3H3","HIST4H4","HJURP","HSP90AA1","HSPA2","HUS1","INCENP","ITGB3BP","KIF18A","KIF20A","KIF23","KIF2A","KIF2B","KIF2C","KNTC1","LIG1","LIN37","LIN52","LIN54","LIN9","LMNA","LMNB1","LOC440577","LOC440917","LOC441488","LOC645084","LOC647654","LOC648152","LOC649620","LOC650621","LOC651610","LOC651763","LOC651921","LOC652826","LOC729964","LOC730418","LOC730594","MAD1L1","MAD2L1","MAPRE1","MAX","MCM10","MCM2","MCM3","MCM4","MCM5","MCM6","MCM7","MCM8","MDM2","MIS12","MIS18A","MIS18BP1","MLF1IP","MNAT1","MYBL2","MYC","NDC80","NDEL1","NEDD1","NEK2","NHP2","NINL","NPM1","NSL1","NUDC","NUF2","NUMA1","NUP107","NUP133","NUP37","NUP43","NUP85","OFD1","OIP5","ORC1","ORC2","ORC3","ORC4","ORC5","ORC6","PAFAH1B1","PCM1","PCNA","PCNT","PKMYT1","PLK1","PLK4","PMF1","POLA1","POLA2","POLD1","POLD2","POLD3","POLD4","POLE","POLE2","POT1","PPP1CC","PPP2CA","PPP2CB","PPP2R1A","PPP2R1B","PPP2R2A","PPP2R3B","PPP2R5A","PPP2R5B","PPP2R5C","PPP2R5D","PPP2R5E","PRIM1","PRIM2","PRKACA","PRKAR2B","PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMA8","PSMB1","PSMB10","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSME1","PSME2","PSME4","PSMF1","PTTG1","RAD1","RAD17","RAD21","RAD9A","RANBP2","RANGAP1","RB1","RBBP4","RBBP4P1","RBBP7","RBL1","RBL2","RCC2","REC8","RFC2","RFC3","RFC4","RFC5","RFWD2","RFWD2P1","RPA1","RPA2","RPA3","RPA4","RPS27","RPS27A","RPS27AP11","RRM2","RSF1","RUVBL1","RUVBL2","SDCCAG8","SEC13","SEH1L","SGOL1","SGOL2","SKA1","SKA2","SKA2L","SKP1","SKP2","SMARCA5","SMC1A","SMC1B","SMC3","SPC24","SPC25","SSNA1","STAG1","STAG2","STAG3","SUN2","SYCP1","SYCP2","SYCP3","SYNE1","SYNE2","TAOK1","TERF1","TERF2","TERF2IP","TERT","TEX12","TFDP1","TINF2","TK2","TP53","TUBA1A","TUBA4A","TUBB","TUBB4A","TUBB4B","TUBBP2","TUBG1","TUBG2","TUBGCP2","TUBGCP3","TUBGCP5","TUBGCP6","TYMS","UBA52","UBE2C","UBE2D1","UBE2E1","UBE2I","WEE1","WRAP53","XPO1","YWHAE","YWHAG","ZW10","ZWILCH","ZWINT")
NSC.genes<-c("NEUROG2","NES","MSI1","SOX2","ASCL1","SOX1","SOX9",
             "FABP7","MSI2","CD133","FGFR4","GLUT1","NEUROD1",
             "RC2","RC1","PAX6","TBR2","EOMES","DCX","BLBP",
             "NCAM","MAP2","TUBB3","FZD9")
TCA.Cycle.genes<-c("ACO2","ADHFE1","BSG","CS","D2HGDH","DLAT","DLD","DLST","FH","IDH1","IDH2","IDH3A","IDH3B","IDH3G","L2HGDH","LDHA","LDHB","LOC283398","LOC646675","LOC646677","LOC650667","LOC650674","LOC650883","LOC651820","MDH2","NNT","OGDH","PDHA1","PDHB","PDHX","PDK1","PDK2","PDK3","PDK4","PDP1","PDP2","PDPR","SDHA","SDHB","SDHC","SDHD","SLC16A1","SLC16A3","SLC16A8","SUCLA2","SUCLA2P1","SUCLG1","SUCLG2")

Neurotransmitter.genes<-c("ABAT","ACHE","ANXA9","BRS3","TSPO","CCKAR","CCKBR","CHAT","CHRM1","CHRM2","CHRM3","CHRNA1","CHRNA2","CHRNA3","CHRNA4","CHRNA5","CHRNA6","CHRNA7","CHRNB1","CHRNB2","CHRNB4","CHRND","CHRNE","CHRNG","COMT","DRD1","DRD2","DRD3","GABRA1","GABRA2","GABRA3","GABRB1","GABRB2","GABRD","GABRE","GABRG1","GABRG2","GABRP","GABRQ","GABRR1","GABRR2","GAD1","GALR1","GALR2","GALR3","GCH1","GCHFR","GLRA1","GLRA2","GLRA3","QRFPR","NPFFR1","MCHR1","PROKR1","PROKR2","NPFFR2","GPR83","GRIA1","GRIN1","GRPR","HCRTR2","HTR1B","HTR2A","HTR3A","HTR3B","MAOA","NMBR","NMUR1","NMUR2","NPY1R","NPY2R","PHOX2A","NPY4R","PRLHR","SLC5A7","SORCS1","SORCS2","SSTR1","SSTR2","SSTR3","SSTR4","TACR1","TACR2","TPH1","B2M","HPRT1","RPL13A","GAPDH","ACTB","HGDC","RTC","RTC","RTC","PPC","PPC","PPC")

Reactome.Apoptosis.genes<-c("ACIN1","ADD1","AKT1","APAF1","APC","APPL1","ARHGAP10","BAD","BAK1","BAX","BBC3","BCAP31","BCL2","BCL2L1","BCL2L11","BID","BIRC2","BMF","BMX","CASP10","CASP3","CASP6","CASP7","CASP8","CASP9","CDH1","CFLAR","CTNNB1","CYCS","DAPK1","DAPK2","DAPK3","DBNL","DCC","DFFA","DFFB","DIABLO","DNM1L","DSG1","DSG2","DSG3","DSP","DYNLL1","DYNLL2","E2F1","FADD","FAS","FASLG","FNTA","GAS2","GSN","GZMB","H1F0","HIST1H1A","HIST1H1B","HIST1H1C","HIST1H1D","HIST1H1E","HMGB1","HMGB2","KPNA1","KPNB1","LMNA","LMNB1","LOC441488","LOC647859","LOC652460","LOC652826","MAGED1","MAPK8","MAPT","MST4","NMT1","OCLN","PAK2","PKP1","PLEC","PMAIP1","PPP3R1","PRKCD","PRKCQ","PSMA1","PSMA2","PSMA3","PSMA4","PSMA5","PSMA6","PSMA7","PSMA8","PSMB1","PSMB10","PSMB2","PSMB3","PSMB4","PSMB5","PSMB6","PSMB7","PSMB8","PSMB9","PSMC1","PSMC2","PSMC3","PSMC4","PSMC5","PSMC6","PSMD1","PSMD10","PSMD11","PSMD12","PSMD13","PSMD14","PSMD2","PSMD3","PSMD4","PSMD5","PSMD6","PSMD7","PSMD8","PSMD9","PSME1","PSME2","PSME4","PSMF1","PTK2","RIPK1","ROCK1","ROCK1P1","RPS27A","RPS27AP11","SATB1","SPTAN1","STK24","TFDP1","TJP1","TJP2","TNF","TNFRSF10B","TNFRSF1A","TNFSF10","TP53","TRADD","TRAF2","UBA52","UNC5A","UNC5B","VIM","XIAP","YWHAB")

GO.Neurogenesis.genes<-c("AGRN","ALS2","AMIGO1","APOE","ARTN","ATP2B2","AZU1","BAI1","BAIAP2","BRSK2","BTG4","CDK5","CDK5R1","CDK6","CIT","CLN5","CNTN4","CYFIP1","DPYSL5","DTX1","EIF2B1","EIF2B2","EIF2B3","EIF2B4","EIF2B5","FARP2","FEZ1","FEZ2","GDNF","GHRL","GLI2","KAL1","KCNIP2","KLK8","KRT2","LAMB1","LDB1","LMX1B","LRRC4C","LST1","MAP1S","MAPT","MDGA1","MDGA2","NF1","NF2","NLGN1","NPTN","NRCAM","NRP1","NRP2","NRTN","NRXN1","NRXN3","NTNG1","NTNG2","OPHN1","OTX2","PARD3","PARD6B","PAX2","PCSK9","PICK1","POU4F1","POU6F2","PPT1","RACGAP1","RND1","ROBO1","ROBO2","RTN1","RTN4","RTN4RL1","RTN4RL2","S100B","SEMA3B","SEMA4F","SERPINF1","SHH","SIAH1","SLIT1","SLIT2","SMARCA1","SOD1","SPON2","TGFB2","THY1","TRAPPC4","UBB","UNC5C","VWC2","YWHAG","YWHAH")

ESC.genes<-c("FOXD3",
             "GATA6",
             "GBX2",
             "NANOG",
             "NR5A2",
             "NR6A1",
             "POU5F1",
             "SOX2",
             "TFCP2L1",
             "UTF1",
             "ZFP42COMMD3",
             "CRABP2",
             "EDNRB",
             "FGF4",
             "FGF5",
             "GABRB3",
             "GAL",
             "GRB7",
             "HCK",
             "IFITM1",
             "IL6ST",
             "KIT",
             "LEFTY1",
             "LEFTY2",
             "LIFR",
             "NODAL",
             "NOG",
             "NUMB",
             "PTEN",
             "SFRP2",
             "TDGF1FGF4",
             "FGF5",
             "GDF3",
             "LEFTY1",
             "LEFTY2",
             "NODAL",
             "TDGF1BRIX1",
             "CD9",
             "DIAPH2",
             "DNMT3B",
             "IFITM2",
             "IGF2BP2",
             "LIN28A",
             "PODXL",
             "REST",
             "SEMA3A",
             "TERTFOXA2",
             "GATA4",
             "PTF1ACDX2",
             "EOMES",
             "GCM1",
             "KRT1AFP",
             "SERPINA1FN1",
             "LAMA1",
             "LAMB1",
             "LAMC1",
             "SOX17T",
             "WT1DES",
             "MYF5",
             "MYOD1HBB",
             "HBZCOL1A1",
             "RUNX2NES",
             "NEUROD1",
             "PAX6CD34",
             "CDH5",
             "FLT1",
             "PECAM1DDX4",
             "SYCP3GCG",
             "IAPP",
             "INS",
             "PAX4",
             "PDX1",
             "SSTOLIG2",
             "TA")

#IUPHAR.genes<-read.delim("IUPHAR_genes.txt",header=F,stringsAsFactors=F)$V1
All_channels<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Ano1","Ano2","Ano3","Ano4","Ano5","Ano6","Ano7","Ano8","Ano9","Ano10","Mip","Aqp1","Aqp2","Aqp3","Aqp4","Aqp5","Aqp6","Aqp7","Aqp8","Aqp9","Aqp10","Aqp11","Aqp12a","Aqp12b","Best1","Best2","Best3","Best4","Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Cftr","Clic1","Clic2","Clic3","Clic4","Clic5","Clic6","Clcnka","Clcnkb","Bsnd","Clcn1","Clcn2","Clcn3","Clcn4","Clcn5","Clcn6","Clcn7","Gja1","Gja3","Gja4","Gja5","Gja6p","Gja8","Gja9","Gja10","Gjb1","Gjb2","Gjb3","Gjb4","Gjb5","Gjb6","Gjb7","Gjc1","Gjc2","Gjc3","Gjd2","Gjd3","Gjd4","Gje1","Itpr1","Itpr2","Itpr3","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Ryr1","Ryr2","Ryr3","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Tpcn1","Tpcn2","Vdac1","Vdac2","Vdac3","Dpp6","Dpp10","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne1b","Kcne2","Kcne3","Kcne4","Kcne5","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4")
Bulk_DE_genes<-c("ENSMUSG00000003309.10","ENSMUSG00000006457.3","ENSMUSG00000015451.13","ENSMUSG00000018470.5","ENSMUSG00000019772.4","ENSMUSG00000019997.8","ENSMUSG00000020178.5","ENSMUSG00000020241.10","ENSMUSG00000020591.10","ENSMUSG00000020878.6","ENSMUSG00000021136.10","ENSMUSG00000022054.8","ENSMUSG00000022123.8","ENSMUSG00000022364.11","ENSMUSG00000022378.10","ENSMUSG00000024085.10","ENSMUSG00000025091.3","ENSMUSG00000026109.11","ENSMUSG00000026765.9","ENSMUSG00000027834.12","ENSMUSG00000028487.15","ENSMUSG00000029813.8","ENSMUSG00000030302.13","ENSMUSG00000030500.5","ENSMUSG00000030583.13","ENSMUSG00000030739.14","ENSMUSG00000031093.11","ENSMUSG00000031410.11","ENSMUSG00000031447.6","ENSMUSG00000031491.9","ENSMUSG00000032076.15","ENSMUSG00000032269.7","ENSMUSG00000032327.11","ENSMUSG00000035131.11","ENSMUSG00000037295.7","ENSMUSG00000037737.6","ENSMUSG00000037852.7","ENSMUSG00000037944.8","ENSMUSG00000038156.12","ENSMUSG00000039391.8","ENSMUSG00000039497.7","ENSMUSG00000040605.6","ENSMUSG00000040759.8","ENSMUSG00000040794.5","ENSMUSG00000042256.4","ENSMUSG00000044748.8","ENSMUSG00000046160.6","ENSMUSG00000047261.9","ENSMUSG00000048138.9","ENSMUSG00000049148.8","ENSMUSG00000054418.6","ENSMUSG00000056158.11","ENSMUSG00000056592.11","ENSMUSG00000059003.9","ENSMUSG00000059327.6","ENSMUSG00000060044.7","ENSMUSG00000061353.8","ENSMUSG00000066705.6","ENSMUSG00000067279.2","ENSMUSG00000074736.7","ENSMUSG00000075014.1","ENSMUSG00000075015.3","ENSMUSG00000075334.2","ENSMUSG00000090919.3","ENSMUSG00000092341.2","ENSMUSG00000094472.1","ENSMUSG00000095280.1","ENSMUSG00000097695.1","ENSMUSG00000097971.3")

Gage_activity_padj_lt_0.01<-c("Arc","Pim1","Ptgs2","Nr4a2","Neat1","Arl4d","Dclk1","Plce1","Baz1a","Fosb","Vgf","Plk2","Atf3","Gpr3","Samd4","Zfp516","Pcdh8","Sgk1","Dpysl5","Siglec1","Lphn3","Midn","Tll1","Phf21b","B230319C09Rik","Egr4","Tsc22d2","Prdm5","Lingo1","Csrnp1","Pcsk1","Ap2b1","Nptx2","Gm25411","Inhba","Baiap2","Bdnf","Spry2","Sertad1","Frmd6","Nr4a3","Nr4a1","Kdm6b","Zdbf2","Blnk","Ddx51","Rnd3","Dmxl1","Calcoco1","Zfp157","Rasl11a","Jdp2","Fbxo33","Homer1","Tacc1","Tmco1","Sik1","Cenpa","Pola2","Fos","Bicd1","Ccrn4l","Rheb","Ankrd33b","Ccdc50","Rfwd2","Nek11","Cstf2","Wdr90","Syne1","Klhdc2","Synpo","Slc25a3","Gm10800","Npas4","Etv3","3110057O12Rik","Cdan1","Hspbap1","Armc9","Xndc1","Gm26870","Gm17024","Dnajc1","Tex29","1700016P03Rik","A730017C20Rik","Cltc","Arsi","Nedd9","Kcnmb4os1","Cep162","Snord16a","Gm12227","Telo2","Sft2d3","Gm3617","Gm26772","Pvrl3","5530601H04Rik","AF357399","Epha5","Tdp1","Rbbp7","Zfand6","Epha7","Ssh2","Clasp1","Prex1","Iqce","Abhd2","Dusp6","B3galt5","Mir673","Snora7a","Bcl7a","Inpp1","Prim2","Mir22hg","Slc6a6","Rasal1","Mir181a-2","Gm28229","Ksr2","Derl2","Hnrnpll","Gm29491","Dusp1","Hltf","Pde4a","Pmepa1","Mir1954","Tmem167","Epha10","Gm12216","Myc","Trip10","Gm17056","Tmem251","Usp45","Ell2","Fmnl1","Amigo3","Dusp5","Dot1l","Rpa1","4931415C17Rik","Glmp","Gm8741","Brip1os","Lmbr1","Ext2","Gm13677","Fam83d","Gm24224","Tet3","Gm13524","Rln3","Gm23202","Ralgps2","Golga5","Exoc1","Capn3","Gm14654","9330159M07Rik","Taok2","Ranbp2","Gm23744","Arid5a","Blvrb","Aldh9a1","Snord12","Mbnl2","Snord65","Zfyve26","Otud3","DQ267101","Scarna3a","Gm27747","Gprc5a","Gpatch2","Oscp1","Csrnp3","Tor1a","Gtf2h2","Shprh","Fhad1os2","Malat1","Tspan33","Gm26767","Blmh","Shc4","1110007C09Rik","Fbxo41","Gem","Gm3076","Tanc1","Mt1","5430405H02Rik","Topbp1","Pld4","Snord49a","Fsd1","Arhgap12","Tas1r3","Fndc3a","Gm16436","Gm6345","Gm10717","Rnpep","Osbpl6","Mir3072","B4galt3","Rgs2","Lmo7","Gadd45g","Gm22358","Adal","Tmem134","Cdkl1","Chst10","Gm6061","Smtn","Gm14239","Gm9917","Armcx1","Gemin8","Tiam1","Fgfr1","Slc26a2","Gm17087","Bai1","Serinc2","Tiam2","Plk3","RP23-459L15.7","Jakmip3","Mcl1","Timm10b","Mthfd2","Ccdc97","Glt8d2","Fam13a","Fam57a","Med14","Bri3","Slc39a3","Hemk1","Gm15201","Pgap1","RP23-382I12.1","Egr3","Kcnj2","Pdk1","Klhl3","Egfl7","Gm28289","Upf1","RP23-14P23.9","Arhgap20","Lrrcc1","Dcaf8","S100a10","Prkab2","Fscn2","Htr1a","Kdm5a","Gm4991","Herc2","Tgfbrap1","Alg6","Acad8","Spag4","Hmga1","Gm26703","Uck1","Klhl5","Farsa","Chrd","Rapgef6","Fam76b","Tmem2","Sds","Tubgcp6","Ncan","Adamts17","Mapk4","Tsen15","Cgref1","Ict1","Paip2b","Clec10a","Ctsz","A830018L16Rik","Gabpb2","Acbd6","Igsf9b","Slc2a1","Map7","Gadd45b","4932438A13Rik","A630023P12Rik","Pkn1","Gm11611","Rbm15","Pbx1","Ift20","Asap2","Ide","Nav1","P4htm","Mapkapk2","Igf2bp2","Brix1","Gm26775","Nub1","Cenpm","Ick","Stxbp4","1500004A13Rik","Orc5","Atg4c","Ptprt","Gm13468","Nek6","Rprl2","Irf3","Gm5601","Mtrr","Mobp","Slc6a17","Gm11794","Dixdc1","Usp30","Stmn4","Gm26530","Fbxl22","Trip12","Rgs7bp","Pam","Flot1","Gm11954","Hspb8","Pms2","Rnu3b4","Rnu3b2","1700102P08Rik","Gpr61","Tarbp2","Npy","Uchl5","Slc6a8","Drc1","Cnot11","Scn3b","Siah2","Gm10801","Scg2","Gm10705","Polg2","Tnfrsf23","Acvr1b","1700019D03Rik","Ifnar1","Kif3c","Cep78","Junb","Cc2d1a","Zfp414","Chgb","Ccnt2","Tyro3","Gm26683","Rell2","Slc16a1","Alpk1","Gm28294","RP24-458F14.4","Skil","Ppm1f","Setdb1","Abca9","Gm20938","Mtfmt","Ankrd32","Dyrk3","Zfp239","Abr","Gprc5c","Lrrc58","D630045J12Rik","RP23-459L15.8","Rrp1b","Stab2","Sbf2","Fbxw7","Gldc","Mysm1","Gpr82","Dirc2","Ptprn","Gls2","Gm15796","Sema3e","Slc6a18","Map3k5","Gm13684","Zfp354c","Zfp831","Ralbp1","Lrrc24","C030039L03Rik","1700064H15Rik","Rasal3","Gm26617","1700001D01Rik","Schip1","Bloc1s2","Nit1","Serpina10","Tmem101","Pak6","A830012C17Rik","Bend4","Gm17229","Kazn","Dcun1d2","Zfp804b","Tsen2")

Gage_FOS_up<-c("Arc","Pim1","Ptgs2","Nr4a2","Neat1","Arl4d","Dclk1","Plce1","Baz1a","Fosb","Vgf","Plk2","Atf3","Gpr3","Samd4","Zfp516","Pcdh8","Sgk1","Dpysl5","Siglec1","Lphn3","Midn","Tll1","Phf21b","B230319C09Rik","Egr4","Tsc22d2","Lingo1","Csrnp1","Pcsk1","Ap2b1","Nptx2","Gm25411","Inhba","Baiap2","Bdnf","Spry2","Sertad1","Frmd6","Nr4a3","Nr4a1","Kdm6b","Zdbf2","Blnk","Ddx51","Rnd3","Dmxl1","Rasl11a","Jdp2","Fbxo33","Homer1","Tacc1","Sik1","Cenpa","Pola2","Fos","Ccrn4l","Rheb","Ankrd33b","Rfwd2","Nek11","Cstf2","Wdr90","Syne1","Synpo","Slc25a3","Gm10800","Npas4","Etv3","Hspbap1","Gm26870","Gm17024","Dnajc1","Tex29","1700016P03Rik","Cltc","Arsi","Nedd9","Kcnmb4os1","Snord16a","Gm12227","Gm3617","Gm26772","Pvrl3","AF357399","Epha5","Rbbp7","Epha7","Ssh2","Prex1","Abhd2","Dusp6","Mir673","Snora7a","Prim2","Mir22hg","Slc6a6","Mir181a-2","Gm28229","Hnrnpll","Gm29491","Dusp1","Pde4a","Pmepa1","Mir1954","Epha10","Myc","Trip10","Gm17056","Tmem251","Ell2","Fmnl1","Amigo3","Dusp5","Dot1l","Rpa1","4931415C17Rik","Gm8741","Brip1os","Gm13677","Fam83d","Gm24224","Tet3","Gm13524","Rln3","Gm23202","Gm14654","9330159M07Rik","Ranbp2","Gm23744","Arid5a","Blvrb","Snord12","Mbnl2","Snord65","Otud3","DQ267101","Scarna3a","Gm27747","Gprc5a","Fhad1os2","Malat1","Gm26767","Blmh","Shc4","1110007C09Rik","Gem","Tanc1","Mt1","Pld4","Snord49a","Tas1r3","Gm16436","Gm6345","Gm10717","Osbpl6","Mir3072","Rgs2","Lmo7","Gadd45g","Gm22358","Gm6061","Smtn","Gm14239","Gm9917","Gemin8","Tiam1","Fgfr1","Gm17087","Bai1","Serinc2","Plk3","RP23-459L15.7","Mcl1","Mthfd2","Med14","Gm15201","Pgap1","RP23-382I12.1","Egr3","Kcnj2","Gm28289","Upf1","S100a10","Prkab2","Fscn2","Kdm5a","Gm4991","Herc2","Spag4","Hmga1","Gm26703","Rapgef6","Tmem2","Sds","Ncan","Mapk4","Cgref1","Clec10a","Ctsz","Igsf9b","Slc2a1","Gadd45b","4932438A13Rik","A630023P12Rik","Gm11611","Rbm15","Asap2","Mapkapk2","Igf2bp2","Brix1","Gm26775","Cenpm","Ptprt","Gm13468","Rprl2","Gm5601","Mobp","Slc6a17","Gm11794","Stmn4","Gm26530","Fbxl22","Trip12","Rgs7bp","Pam","Hspb8","Rnu3b4","Rnu3b2","1700102P08Rik","Npy","Uchl5","Slc6a8","Cnot11","Siah2","Gm10801","Scg2","Polg2","Tnfrsf23","Junb","Chgb","Ccnt2","Tyro3","Gm26683","Rell2","Slc16a1","Alpk1","Gm28294","RP24-458F14.4","Skil","Abca9","Gm20938","Dyrk3","Abr","Lrrc58","RP23-459L15.8","Rrp1b","Stab2","Fbxw7","Gldc","Gpr82","Dirc2","Ptprn","Gls2","Gm15796","Sema3e","Slc6a18","Map3k5","Gm13684","1700064H15Rik","Rasal3","Gm26617","1700001D01Rik","Schip1","Serpina10","Pak6","A830012C17Rik","Bend4","Gm17229","Zfp804b")
Gage_FOS_dn<-c("Prdm5","Calcoco1","Zfp157","Tmco1","Bicd1","Ccdc50","Klhdc2","3110057O12Rik","Cdan1","Armc9","Xndc1","A730017C20Rik","Cep162","Telo2","Sft2d3","5530601H04Rik","Tdp1","Zfand6","Clasp1","Iqce","B3galt5","Bcl7a","Inpp1","Rasal1","Ksr2","Derl2","Hltf","Tmem167","Gm12216","Usp45","Glmp","Lmbr1","Ext2","Ralgps2","Golga5","Exoc1","Capn3","Taok2","Aldh9a1","Zfyve26","Gpatch2","Oscp1","Csrnp3","Tor1a","Gtf2h2","Shprh","Tspan33","Fbxo41","Gm3076","5430405H02Rik","Topbp1","Fsd1","Arhgap12","Fndc3a","Rnpep","B4galt3","Adal","Tmem134","Cdkl1","Chst10","Armcx1","Slc26a2","Tiam2","Jakmip3","Timm10b","Ccdc97","Glt8d2","Fam13a","Fam57a","Bri3","Slc39a3","Hemk1","Pdk1","Klhl3","Egfl7","RP23-14P23.9","Arhgap20","Lrrcc1","Dcaf8","Htr1a","Tgfbrap1","Alg6","Acad8","Uck1","Klhl5","Farsa","Chrd","Fam76b","Tubgcp6","Adamts17","Tsen15","Ict1","Paip2b","A830018L16Rik","Gabpb2","Acbd6","Map7","Pkn1","Pbx1","Ift20","Ide","Nav1","P4htm","Nub1","Ick","Stxbp4","1500004A13Rik","Orc5","Atg4c","Nek6","Irf3","Mtrr","Dixdc1","Usp30","Flot1","Gm11954","Pms2","Gpr61","Tarbp2","Drc1","Scn3b","Gm10705","Acvr1b","1700019D03Rik","Ifnar1","Kif3c","Cep78","Cc2d1a","Zfp414","Ppm1f","Setdb1","Mtfmt","Ankrd32","Zfp239","Gprc5c","D630045J12Rik","Sbf2","Mysm1","Zfp354c","Zfp831","Ralbp1","Lrrc24","C030039L03Rik","Bloc1s2","Nit1","Tmem101","Kazn","Dcun1d2","Tsen2")

GluR_genes<-c("Gria1","Gria2","Gria3","Gria4","Grid1","Grid2","Grik1","Grik2","Grik3","Grik4","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b","Grm1","Grm2","Grm3","Grm4","Grm5","Grm6","Grm7","Grm8","Grina")

Ligand_gated_ion_channels<-c("Htr3a","Htr3b","Htr3c","Htr3d","Htr3e","Cftr","Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10","Chrnb1","Chrnb2","Chrnb3","Chrnb4","Chrnd","Chrne","Chrng","Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6","Gabrb1","Gabrb2","Gabrb3","Gabrd","Gabre","Gabrg1","Gabrg2","Gabrg3","Gabrp","Gabrq","Gabrr1","Gabrr2","Gabrr3","Gria1","Gria2","Gria3","Gria4","Grid1","Grid2","Grik1","Grik2","Grik3","Grik4","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b","Glra1","Glra2","Glra3","Glra4","Glrb","Itpr1","Itpr2","Itpr3","P2rx1","P2rx2","P2rx3","P2rx4","P2rx5","P2rx6","P2rx7","Ryr1","Ryr2","Ryr3","Zacn")

Synaptic_transmission_genes<-c("2610042L04Rik","9530002B09Rik","Abat","Abhd6","Adcyap1","Adgrl1","Adipoq","Adnp","Adora1","Adora2a","Adra1a","Adrb1","Adrb2","Agrn","Agt","Als2","Anapc2","Apba1","Apba2","Apba3","Arc","Arf1","Arrb2","Asic1","Atad1","Atp2a2","Atp2b2","Atxn1","Baiap2","Bche","Bdnf","Braf","Brsk1","Btbd9","Cacna1a","Cacna1b","Cacna1c","Cacna1g","Cacna2d2","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng2","Cacng3","Cacng4","Cacng5","Cacng7","Cacng8","Cadps","Cadps2","Calb1","Camk2a","Camk2b","Camk4","Car2","Car7","Cartpt","Cav2","Cckbr","Ccl2","Cd24a","Cd38","Cdc20","Cdh8","Cdk5","Cel","Celf4","Celsr1","Chat","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10","Chrnb1","Chrnb2","Chrnb3","Chrnb4","Chrnd","Chrne","Chrng","Clcn3","Clstn1","Clstn2","Clstn3","Cnih2","Cnih3","Cnr1","Cnr2","Cntn2","Cntnap4","Cpeb1","Cpeb3","Cplx1","Cplx2","Cplx3","Cplx4","Creb1","Crh","Crhbp","Crhr1","Crhr2","Crtc1","Cspg5","Ctnnb1","Ctnnd2","Cux2","Cx3cl1","Cx3cr1","D730048I06Rik","Dapk1","Dbi","Dbn1","Dgki","Dkk1","Dlg1","Dlg2","Dlg3","Dlg4","Dlgap1","Dlgap2","Dmpk","Dnm1","Doc2a","Doc2b","Doc2g","Drd1","Drd2","Drd3","Drd4","Drd5","Dtnbp1","Dvl1","Edn1",
  "Edn3","Egfr","Egr1","Egr3","Eif2ak4","Eif4a3","Eif4ebp2","Ephb2","Erbb4","Etv5","Exoc4","Fgf12","Fgf14","Flot1","Fmr1","Gabbr1","Gabra1","Gabra3","Gabra5","Gabra6","Gabrb2","Gabrb3","Gabrd","Gabrg1","Gabrg2","Gdnf","Gfap","Ghrl","Gip","Gipc1","Gjc1","Gjd2","Glra1","Glra2","Glra3","Glra4","Glrb","Gls","Glul","Gm12886","Gm12887","Gm12888","Gnai1","Gnai2","Gpm6b","Gpr1","Gpr21","Gpr52","Gpr149","Gria1","Gria2","Gria4","Grid2","Grid2ip","Grik1","Grik2","Grik3","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grm1","Grm2","Grm3","Grm4","Grm5","Grm6","Grm7","Grm8","Gsg1l","Gsk3b","Hap1","Hcrt","Hcrtr1","Hcrtr2","Hnrnpk","Hras","Hrh1","Hrh2","Hrh3","Hrh4","Htr1b","Htr1d","Htr2a","Htr2c","Htr4","Htr6","Htr7","Htt","Ica1","Ifng","Igsf9b","Iqsec2","Itpka","Itpr3","Jph3","Jph4","Kalrn","Kcnb1","Kcnc3","Kcnc4","Kcnj10","Kcnma1","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4","Kcnn2","Kcnq4","Kdr","Kif1b","Kif5b","Kiss1","Kiss1r","Kit","Klhl24","Kmt2a","Kras","Lama2","Lgi1","Lin7a","Lin7b","Lin7c","Lrp6","Lrp8","Lrrk2","Lrrtm1","Lrrtm2","Lynx1","Lypd1","Lzts1","Mapk1","Mapk8ip2","Mapt","Mecp2","Mef2c","Met","Mgll","Mink1","Mmp9","Mtmr2","Musk","Mylk2","Myo5a","Myo6","Napa","Napb","Nat8l","Ncam1","Ncdn","Ncs1","Neto1","Neto2","Neurl1a","Neurod2","Nf1","Nfatc4","Ngf","Ngfr","Nlgn1","Nlgn2","Nlgn3","Nlgn4l","Nmu","Nos1","Npas4","Npbwr1","Npff","Npffr1","Npffr2","Nps","Nptn","Npy2r","Npy5r","Nr2e1",
  "Nr3c1","Nrgn","Nrxn1","Nrxn2","Nrxn3","Nsmf","Ntf3","Ntrk1","Ntrk2","Ntsr1","Omp","Ophn1","Oprd1","Oprk1","Oprl1","Oprm1","Otof","Oxt","Oxtr","P2rx1","P2rx2","P2rx3","P2rx4","P2rx7","Pafah1b1","Paip2","Pak1","Park2","Park7","Pate4","Pcdh8","Pcdh17","Pcdhb16","Pclo","Pdyn","Pdzd11","Pebp1","Penk","Pfn2","Pick1","Pink1","Pla2g6","Plat","Plcb3","Plcl1","Plcl2","Plk2","Pmch","Pnkd","Pnoc","Ppargc1a","Ppfia3","Ppp1r9a","Ppp1r9b","Ppp3ca","Ppt1","Prkaca","Prkce","Prkcg","Prkcz","Psen1","Psen2","Pten","Ptgdr2","Ptger4","Ptgs2","Ptk2","Ptk2b","Ptpn5","Ptprn2","Rab3a","Rab3b","Rab3gap1","Rab5a","Rab8a","Rab11a","Rac1","Rac3","Rap1a","Rap1b","Rapgef2","Rapgef4","Rapsn","Rara","Rasd2","Rasgrf1","Rasgrf2","Reln","Retn","Rgs14","Ric3","Rims1","Rims2","Rims3","Rims4","Rin1","Rph3a","S1pr2","S100b","Scrib","Sdcbp","Serpine2","Sez6","Shank1","Shank2","Shank3","Shc3","Shisa6","Shisa7","Shisa8","Shisa9","Shisa9","Sipa1l1","Slc1a3","Slc5a7","Slc6a1","Slc6a3","Slc6a4","Slc6a5","Slc6a9","Slc6a11","Slc6a12","Slc6a13","Slc6a20a","Slc6a20b","Slc8a2","Slc8a3","Slc12a4","Slc12a5","Slc12a6","Slc12a7","Slc17a7","Slc18a1","Slc24a2","Slc29a1","Slc30a1","Slitrk5","Snap23","Snap25","Snap29","Snap47","Snap91",
  "Snapin","Snca","Sncaip","Sncb","Sncg","Sorcs3","Spg11","Sphk1","Sptbn2","Srf","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","Stac3","Star","Stau1","Stau2","Ston2","Stx1a","Stx1b","Stx2","Stx3","Stx4a","Stx11","Stx19","Stxbp1","Sv2a","Sv2b","Sv2c","Syde1","Syn1","Syn2","Syn3","Syngap1","Syngr1","Synj1","Syp","Syt1","Syt2","Syt3","Syt4","Syt5","Syt6","Syt7","Syt8","Syt9","Syt10","Syt11","Syt12","Syt13",
  "Syt15","Syt17","Sytl1","Sytl3","Sytl4","Sytl5","Tac1","Tacr1","Tacr2","Th","Tmod2","Tnf","Tnr","Tor1a","Tpgs1","Trim9","Ucn","Unc13a","Unc13b","Unc13c","Usp14","Usp46","Uts2","Vamp2","Vdac1","Vdac3","Vgf","Wnt7a","Xbp1","Ywhag","Zmynd8")

Amine_receptors<-c("Htr1a","Htr1b","Htr1d","Htr1e","Htr1f","Htr2a","Htr2b","Htr2c","Htr4","Htr5a","Htr5bp","Htr6","Htr7","Adra1a","Adra1b","Adra1d","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Adrb3","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","Drd1","Drd2","Drd3","Drd4","Drd5","Hrh1","Hrh2","Hrh3","Hrh4","Taar1","Taar2","Taar3","Taar4p","Taar5","Taar6","Taar7p","Taar8","Taar9")

Peptide_receptors<-c("Aplnr","Ghsr","Kiss1r","Mlnr","Prlhr","Qrfpr","Trhr","Agtr1","Agtr2","Avpr1a","Avpr1b","Avpr2","Oxtr","Bdkrb1","Bdkrb2","Cckar","Cckbr","Ednra","Ednrb","Galr1","Galr2","Galr3","Gnrhr","Gnrhr2","Hcrtr1","Hcrtr2","Mchr1","Mchr2","Mc1r","Mc2r","Mc3r","Mc4r","Mc5r","Nmur1","Nmur2","Npbwr1","Npbwr2","Npffr1","Npffr2","Npsr1","Npy1r","Npy2r","Npy4r","Npy5r","Npy6r","Ntsr1","Ntsr2","Oprd1","Oprk1","Oprl1","Oprm1","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","Tacr1","Tacr2","Tacr3")

VIP_calcitonin<-c("Adcyap1r1","Vipr1","Vipr2","Calcr","Calcrl","Crhr1","Crhr2","Casr","Gprc6a")

Other_receptors<-c(Amine_receptors,Peptide_receptors,VIP_calcitonin,"Pgrmc1","Pgrmc2","Nenf","Cyb5d2","Glp2r")

GPCRs<-c("Htr1a","Htr1b","Htr1d","Htr1e","Htr1f","Htr2a","Htr2b","Htr2c","Htr4","Htr5a","Htr5bp","Htr6","Htr7","Gpr107","Gpr108","Gpr137","Gpr143","Gpr157","Tpra1","Adora1","Adora2a","Adora2b","Adora3","Adgra1","Adgra2","Adgra3","Adgrb1","Adgrb2","Adgrb3","Celsr1","Celsr2","Celsr3","Adgrd1","Adgrd2","Adgre1","Adgre2","Adgre3","Adgre4p","Adgre5","Adgrf1","Adgrf2","Adgrf3","Adgrf4","Adgrf5","Adgrg1","Adgrg2","Adgrg3","Adgrg4","Adgrg5","Adgrg6","Adgrg7","Adgrl1","Adgrl2","Adgrl3","Adgrl4","Adgrv1","Adra1a","Adra1b","Adra1d","Adra2a","Adra2b","Adra2c","Adrb1","Adrb2","Adrb3","Agtr1","Agtr2","Avpr1a","Avpr1b","Avpr2","Oxtr","Ackr1","Ackr2","Ackr3","Ackr4","Ccrl2","Pitpnm3","Bdkrb1","Bdkrb2","Ccr1","Ccr10","Ccr2","Ccr3","Ccr4","Ccr5","Ccr6","Ccr7","Ccr8","Ccr9","Cx3cr1","Cxcr1","Cxcr2","Cxcr3"
         ,"Cxcr4","Cxcr5","Cxcr6","Calcr","Calcrl","Casr","Gprc6a","Cnr1","Cnr2","Cmklr1","Cckar","Cckbr","Chrm1","Chrm2","Chrm3","Chrm4","Chrm5","C3ar1","C5ar1","C5ar2","Crhr1","Crhr2","Drd1","Drd2","Drd3","Drd4","Drd5","Ednra","Ednrb","F2r","F2rl1","F2rl2","F2rl3","Fpr1","Fpr2","Fpr3","Ffar1","Ffar2","Ffar3","Ffar4","Gpbar1","Gper1","Gpr1","Gpr101","Gpr119","Gpr12","Gpr132","Gpr135","Gpr139","Gpr141","Gpr142","Gpr146","Gpr148","Gpr149","Gpr15","Gpr150","Gpr151","Gpr152","Gpr153","Gpr160","Gpr161","Gpr162","Gpr17","Gpr171","Gpr173","Gpr174","Gpr176","Gpr18","Gpr182","Gpr183","Gpr19","Gpr20","Gpr21","Gpr22","Gpr25","Gpr26","Gpr27","Gpr3","Gpr31","Gpr32","Gpr33","Gpr34","Gpr35","Gpr37","Gpr37l1","Gpr39","Gpr4","Gpr42","Gpr45","Gpr50","Gpr52","Gpr55","Gpr6","Gpr61","Gpr62","Gpr63","Gpr65","Gpr68","Gpr75","Gpr78","Gpr79","Gpr82","Gpr83","Gpr84","Gpr85","Gpr87","Gpr88","Lgr4","Lgr5","Lgr6","Mas1","Mas1l","Mrgprd","Mrgpre","Mrgprf","Mrgprg","Mrgprx1","Mrgprx2","Mrgprx3","Mrgprx4","Gpr156","Gpr158","Gpr179","Gprc5a","Gprc5b","Gprc5c","Gprc5d","Fzd1","Fzd10","Smo","Fzd2","Fzd3","Fzd4","Fzd5","Fzd6","Fzd7","Fzd8","Fzd9","Galr1","Galr2","Galr3","Gabbr1","Gabbr2","Gcgr","Ghrhr","Gipr","Glp1r","Glp2r","Sctr","Grm1","Grm2","Grm3","Grm4","Grm5","Grm6","Grm7","Grm8","Fshr","Lhcgr","Tshr","Gnrhr","Gnrhr2","Hrh1","Hrh2","Hrh3","Hrh4","Hcar1","Hcar2","Hcar3","Hcrtr1","Hcrtr2","Cysltr1","Cysltr2","Fpr2","Ltb4r","Ltb4r2","Oxer1","Lpar1","Lpar2",
         "Lpar3","Lpar4","Lpar5","Lpar6","Mchr1","Mchr2","Mc1r","Mc2r","Mc3r","Mc4r","Mc5r","Mtnr1a","Mtnr1b","Nmur1","Nmur2","Npbwr1","Npbwr2","Npffr1","Npffr2","Npsr1","Npy1r","Npy2r","Npy4r","Npy5r","Npy6r","Ntsr1","Ntsr2","Or1a1","Or1a2","Or1aa1p","Or1ab1p","Or1ac1p","Or1b1","Or1c1","Or1d2","Or1d3p","Or1d4","Or1d5","Or1e1","Or1e2","Or1e3","Or1f1","Or1f12","Or1f2p","Or1g1","Or1h1p","Or1i1","Or1j1","Or1j2","Or1j4","Or1k1","Or1l1","Or1l3","Or1l4","Or1l6","Or1l8","Or1m1","Or1m4p","Or1n1","Or1n2","Or1p1","Or1q1","Or1r1p","Or1s1","Or1s2","Or1x1p","Or1x5p","Or10a2","Or10a3","Or10a4","Or10a5","Or10a6","Or10a7","Or10aa1p","Or10ab1p","Or10ac1","Or10ad1","Or10ae1p","Or10ae3p","Or10af1p","Or10ag1","Or10ah1p","Or10ak1p","Or10b1p","Or10c1","Or10d1p","Or10d3","Or10d4p","Or10d5p","Or10g1p","Or10g2","Or10g3","Or10g4","Or10g5p","Or10g6","Or10g7","Or10g8","Or10g9","Or10h1","Or10h2","Or10h3","Or10h4","Or10h5","Or10j1","Or10j2p","Or10j3","Or10j4","Or10j5","Or10j6p","Or10j7p","Or10j8p","Or10j9p","Or10k1","Or10k2","Or10n1p","Or10p1","Or10q1","Or10q2p","Or10r1p","Or10r2","Or10r3p","Or10s1","Or10t1p","Or10t2","Or10u1p","Or10v1","Or10v2p","Or10v3p","Or10v7p","Or10w1","Or10x1","Or10y1p","Or10z1","Or11a1","Or11g1p","Or11g2","Or11h1","Or11h12","Or11h13p","Or11h2","Or11h3p","Or11h4","Or11h5p","Or11h6","Or11h7","Or11i1p","Or11j1p","Or11j2p","Or11j5p","Or11k1p","Or11k2p","Or11l1","Or11m1p","Or11n1p","Or11p1p","Or11q1p","Or12d1","Or12d2","Or12d3","Or13a1"
         ,"Or13c1p","Or13c2","Or13c3","Or13c4","Or13c5","Or13c6p","Or13c7","Or13c8","Or13c9","Or13d1","Or13d2p","Or13d3p","Or13e1p","Or13f1","Or13g1","Or13h1","Or13i1p","Or13j1","Or13k1p","Or13z1p","Or13z2p","Or13z3p","Or14a16","Or14a2","Or14c36","Or14i1","Or14j1","Or14k1","Or14l1p","Or2a1","Or2a12","Or2a13p","Or2a14","Or2a15p","Or2a2","Or2a20p","Or2a25","Or2a3p","Or2a4","Or2a41p","Or2a42","Or2a5","Or2a7","Or2a9p","Or2ad1p","Or2ae1","Or2af1p","Or2ag1","Or2ag2","Or2ah1p","Or2ai1p","Or2aj1","Or2ak2","Or2al1p","Or2am1p","Or2ao1p","Or2ap1","Or2aq1p","Or2as1p","Or2as2p","Or2at1p","Or2at2p","Or2at4","Or2b11","Or2b2","Or2b3","Or2b4p","Or2b6","Or2b7p","Or2b8p","Or2bh1p","Or2c1","Or2c3","Or2d2","Or2d3","Or2e1p","Or2f1","Or2f2","Or2g1p","Or2g2","Or2g3","Or2g6","Or2h1","Or2h2","Or2h4p","Or2h5p","Or2i1p","Or2j1","Or2j2","Or2j3","Or2j4p","Or2k2","Or2l13","Or2l1p","Or2l2","Or2l3","Or2l5","Or2l6p","Or2l8","Or2l9p","Or2m1p","Or2m2","Or2m3","Or2m4","Or2m5","Or2m7","Or2n1p","Or2p1p","Or2q1p","Or2r1p","Or2s1p","Or2s2","Or2t1","Or2t10","Or2t11","Or2t12","Or2t2","Or2t27","Or2t29","Or2t3","Or2t32p","Or2t33","Or2t34","Or2t35","Or2t4","Or2t5","Or2t6","Or2t7","Or2t8","Or2u1p","Or2u2p","Or2v1","Or2v2","Or2w1","Or2w2p","Or2w3","Or2w4p","Or2w5","Or2w6p","Or2x1p","Or2y1","Or2z1","Or3a1","Or3a2","Or3a3","Or3a4p","Or3b1p","Or3d1p","Or4a10p","Or4a11p","Or4a12p","Or4a13p","Or4a14p","Or4a15","Or4a16","Or4a17p","Or4a18p","Or4a19p","Or4a1p","Or4a21p","Or4a2p","Or4a3p","Or4a40p","Or4a41p","Or4a42p","Or4a43p","Or4a44p","Or4a45p","Or4a46p","Or4a47","Or4a48p","Or4a49p","Or4a4p","Or4a5","Or4a50p","Or4a6p","Or4a7p","Or4a8",
         "Or4a9p","Or4b1","Or4b2p","Or4c10p","Or4c11","Or4c12","Or4c13","Or4c14p","Or4c15","Or4c16","Or4c1p","Or4c2p","Or4c3","Or4c45","Or4c46","Or4c48p","Or4c49p","Or4c4p","Or4c5","Or4c50p","Or4c6","Or4c7p","Or4c9p","Or4d1","Or4d10","Or4d11","Or4d12p","Or4d2","Or4d5","Or4d6","Or4d7p","Or4d8p","Or4d9","Or4e1","Or4e2","Or4f13p","Or4f14p","Or4f15","Or4f16","Or4f17","Or4f1p","Or4f21","Or4f28p","Or4f29","Or4f2p","Or4f3","Or4f4","Or4f5","Or4f6","Or4f7p","Or4f8p","Or4g11p","Or4g1p","Or4g2p","Or4g3p","Or4g4p","Or4g6p","Or4h12p","Or4h6p","Or4k1","Or4k11p","Or4k12p","Or4k13","Or4k14","Or4k15","Or4k16p","Or4k17","Or4k2","Or4k3","Or4k4p","Or4k5","Or4k6p","Or4k7p","Or4k8p","Or4l1","Or4m1","Or4m2","Or4n1p","Or4n2","Or4n3p","Or4n4","Or4n5","Or4p1p","Or4p4","Or4q1p","Or4q2","Or4q3","Or4r1p","Or4r2p","Or4r3p","Or4s1","Or4s2","Or4t1p","Or4u1p","Or4v1p","Or4w1p","Or4x1","Or4x2","Or4x7p","Or5a1","Or5a2","Or5ac1","Or5ac2","Or5ac4p","Or5ah1p","Or5ak1p","Or5ak2","Or5ak3p","Or5ak4p","Or5al1","Or5al2p","Or5am1p","Or5an1","Or5an2p","Or5ao1p","Or5ap1p","Or5ap2","Or5aq1p","Or5ar1","Or5as1","Or5au1","Or5aw1p","Or5az1p","Or5b10p","Or5b12","Or5b15p","Or5b17","Or5b19p","Or5b1p","Or5b2","Or5b21","Or5b3","Or5ba1p","Or5bb1p","Or5bc1p","Or5bd1p","Or5be1p","Or5bh1p","Or5bj1p","Or5bk1p","Or5bl1p","Or5bm1p"
         ,"Or5bn1p","Or5bn2p","Or5bp1p","Or5bq1p","Or5br1p","Or5bs1p","Or5bt1p","Or5c1","Or5d13","Or5d14","Or5d15p","Or5d16","Or5d17p","Or5d18","Or5d2p","Or5d3p","Or5e1p","Or5f1","Or5f2p","Or5g1p","Or5g3","Or5g4p","Or5g5p","Or5h1","Or5h14","Or5h15","Or5h2","Or5h3p","Or5h4p","Or5h5p","Or5h6","Or5h7p","Or5h8","Or5i1","Or5j1p","Or5j2","Or5j7p","Or5k1","Or5k2","Or5k3","Or5k4","Or5l1","Or5l2","Or5m1","Or5m10","Or5m11","Or5m12p","Or5m13p","Or5m14p","Or5m2p","Or5m3","Or5m4p","Or5m5p","Or5m6p","Or5m7p","Or5m8","Or5m9","Or5p1p","Or5p2","Or5p3","Or5p4p","Or5r1","Or5s1p","Or5t1","Or5t2","Or5t3","Or5v1","Or5w1p","Or5w2","Or51a10p","Or51a1p","Or51a2","Or51a3p","Or51a4","Or51a5p","Or51a6p","Or51a7","Or51a8p","Or51a9p","Or51ab1p","Or51b2","Or51b3p","Or51b4","Or51b5","Or51b6","Or51b8p","Or51c1p","Or51c4p","Or51d1","Or51e1","Or51e2","Or51f1","Or51f2","Or51f3p","Or51f4p","Or51f5p","Or51g1","Or51g2","Or51h1","Or51h2p","Or51i1","Or51i2","Or51j1","Or51k1p","Or51l1","Or51m1","Or51n1p","Or51p1p","Or51q1","Or51r1p","Or51s1","Or51t1","Or51v1","Or52a1","Or52a4p","Or52a5","Or52b1p","Or52b2","Or52b3p","Or52b4","Or52b5p","Or52b6","Or52d1","Or52e1","Or52e2","Or52e3p","Or52e4","Or52e5","Or52e6","Or52e7p","Or52e8","Or52h1","Or52h2p","Or52i1","Or52i2","Or52j1p","Or52j2p","Or52j3","Or52k1","Or52k2","Or52k3p","Or52l1","Or52l2p","Or52m1","Or52m2p","Or52n1","Or52n2","Or52n3p","Or52n4","Or52n5","Or52p1p","Or52p2p","Or52q1p","Or52r1","Or52s1p","Or52t1p","Or52u1p",
         "Or52v1p","Or52w1","Or52x1p","Or52y1p","Or52z1","Or55b1p","Or56a1","Or56a3","Or56a4","Or56a5","Or56a7p","Or56b1","Or56b2p","Or56b3p","Or56b4","Or6a2","Or6b1","Or6b2","Or6b3","Or6c1","Or6c2","Or6c3","Or6c4","Or6c5p","Or6c6","Or6c64p","Or6c65","Or6c66p","Or6c68","Or6c69p","Or6c70","Or6c71p","Or6c72p","Or6c73p","Or6c74","Or6c75","Or6c76","Or6c7p","Or6d1p","Or6e1p","Or6f1","Or6j1","Or6k1p","Or6k2","Or6k3","Or6k4p","Or6k5p","Or6k6","Or6l1p","Or6l2p","Or6m1","Or6m2p","Or6m3p","Or6n1","Or6n2","Or6p1","Or6q1","Or6r1p","Or6r2p","Or6s1","Or6t1","Or6u2p","Or6v1","Or6w1p","Or6x1","Or6y1","Or7a10","Or7a11p","Or7a15p","Or7a17","Or7a18p","Or7a19p","Or7a1p","Or7a2p","Or7a3p","Or7a5","Or7a8p","Or7c1","Or7c2","Or7d11p","Or7d1p","Or7d2","Or7d4","Or7e100p","Or7e101p","Or7e102p","Or7e104p","Or7e105p","Or7e106p","Or7e108p","Or7e109p","Or7e10p","Or7e110p","Or7e111p","Or7e115p","Or7e116p","Or7e117p","Or7e11p","Or7e121p","Or7e122p","Or7e125p","Or7e126p","Or7e128p","Or7e129p","Or7e12p","Or7e130p","Or7e136p","Or7e13p","Or7e140p","Or7e145p","Or7e148p","Or7e149p","Or7e14p","Or7e154p","Or7e155p","Or7e156p","Or7e157p","Or7e158p","Or7e159p","Or7e15p","Or7e160p","Or7e161p","Or7e162p","Or7e163p","Or7e16p","Or7e18p","Or7e19p","Or7e1p","Or7e21p","Or7e22p","Or7e23p","Or7e24","Or7e25p","Or7e26p","Or7e28p","Or7e29p","Or7e2p","Or7e31p","Or7e33p","Or7e35p","Or7e36p","Or7e37p","Or7e38p","Or7e39p","Or7e41p","Or7e43p","Or7e46p","Or7e47p","Or7e4p","Or7e53p","Or7e55p","Or7e59p","Or7e5p","Or7e62p","Or7e66p","Or7e7p","Or7e83p","Or7e84p","Or7e85p","Or7e86p","Or7e87p","Or7e89p","Or7e8p","Or7e90p","Or7e91p","Or7e93p","Or7e94p","Or7e96p","Or7e97p","Or7e99p","Or7g1","Or7g15p","Or7g2","Or7g3","Or7h1p","Or7h2p","Or7k1p","Or7l1p","Or7m1p","Or8a1","Or8a2p","Or8a3p","Or8b10p","Or8b12","Or8b1p","Or8b2","Or8b3","Or8b4","Or8b5p","Or8b6p","Or8b7p","Or8b8","Or8b9p","Or8c1p","Or8d1","Or8d2","Or8d4","Or8f1p","Or8g1","Or8g2","Or8g3p","Or8g5","Or8g7p","Or8h1","Or8h2","Or8h3","Or8i1p","Or8i2","Or8i4p","Or8j1","Or8j2","Or8j3","Or8k1","Or8k2p","Or8k3","Or8k4p","Or8k5","Or8l1p","Or8q1p","Or8r1p","Or8s1","Or8s21p","Or8t1p","Or8u1","Or8u8","Or8u9","Or8v1p","Or8x1p","Or9a1p","Or9a2","Or9a3p","Or9a4","Or9g1","Or9g2p","Or9g3p","Or9g4","Or9g9","Or9h1p","Or9i1","Or9i2p","Or9i3p"
         ,"Or9k1p","Or9k2","Or9l1p","Or9m1p","Or9n1p","Or9p1p","Or9q1","Or9q2","Or9r1p","Or9s24p","Oprd1","Oprk1","Oprl1","Oprm1","Opn1lw","Opn1mw","Opn1mw2","Opn1mw3","Opn1sw","Rho","Opn3","Opn4","Opn5","Rgr","Rrh","Oxgr1","Pth1r","Pth2r","Aplnr","Ghsr","Kiss1r","Mlnr","Prlhr","Qrfpr","Trhr","Ptafr","Prokr1","Prokr2","Ptgdr","Ptgdr2","Ptger1","Ptger2","Ptger3","Ptger4","Ptgfr","Ptgir","Tbxa2r","P2ry1","P2ry10","P2ry11","P2ry12","P2ry13","P2ry14","P2ry2","P2ry4","P2ry6","P2ry8","Rxfp1","Rxfp2","Rxfp3","Rxfp4","Sstr1","Sstr2","Sstr3","Sstr4","Sstr5","S1pr1","S1pr2","S1pr3","S1pr4","S1pr5","Sucnr1","Tacr1","Tacr2","Tacr3","Tas1r1","Tas1r2","Tas1r3","Tas2r1","Tas2r10","Tas2r12p","Tas2r13","Tas2r14","Tas2r15p","Tas2r16","Tas2r18p","Tas2r19","Tas2r20","Tas2r22","Tas2r2p","Tas2r3","Tas2r30","Tas2r31","Tas2r33","Tas2r36","Tas2r37","Tas2r38","Tas2r39","Tas2r4","Tas2r40","Tas2r41","Tas2r42","Tas2r43","Tas2r45","Tas2r46","Tas2r5","Tas2r50","Tas2r60","Tas2r62p","Tas2r63p","Tas2r64p","Tas2r67p","Tas2r68p","Tas2r6p","Tas2r7","Tas2r8","Tas2r9","Taar1","Taar2","Taar3","Taar4p","Taar5","Taar6","Taar7p","Taar8","Taar9","Adcyap1r1","Vipr1","Vipr2","Vn1r1","Vn1r100p","Vn1r101p","Vn1r102p","Vn1r103p","Vn1r104p","Vn1r105p","Vn1r106p","Vn1r107p","Vn1r108p","Vn1r109p","Vn1r10p","Vn1r110p","Vn1r111p","Vn1r112p","Vn1r11p","Vn1r12p","Vn1r13p","Vn1r14p","Vn1r15p","Vn1r16p","Vn1r17p","Vn1r18p","Vn1r19p","Vn1r2","Vn1r20p","Vn1r21p","Vn1r22p","Vn1r23p","Vn1r24p","Vn1r25p","Vn1r26p","Vn1r27p","Vn1r28p","Vn1r29p","Vn1r3","Vn1r30p","Vn1r31p","Vn1r32p","Vn1r33p","Vn1r34p","Vn1r35p","Vn1r36p","Vn1r37p","Vn1r38p","Vn1r39p","Vn1r4","Vn1r40p","Vn1r41p","Vn1r42p","Vn1r43p","Vn1r44p","Vn1r45p","Vn1r46p","Vn1r47p","Vn1r48p","Vn1r49p","Vn1r5","Vn1r51p","Vn1r52p","Vn1r53p","Vn1r54p"
         ,"Vn1r55p","Vn1r56p","Vn1r57p","Vn1r58p","Vn1r59p","Vn1r60p","Vn1r61p","Vn1r62p","Vn1r63p","Vn1r64p","Vn1r65p","Vn1r66p","Vn1r67p","Vn1r68p","Vn1r69p","Vn1r6p","Vn1r70p","Vn1r71p","Vn1r72p","Vn1r73p","Vn1r74p","Vn1r75p","Vn1r76p","Vn1r77p","Vn1r78p","Vn1r79p","Vn1r7p","Vn1r80p","Vn1r81p","Vn1r82p","Vn1r83p","Vn1r84p","Vn1r85p","Vn1r86p","Vn1r87p","Vn1r88p","Vn1r89p","Vn1r8p","Vn1r90p","Vn1r91p","Vn1r92p","Vn1r93p","Vn1r94p","Vn1r95p","Vn1r96p","Vn1r97p","Vn1r98p","Vn1r99p","Vn1r9p","Vn2r10p","Vn2r11p","Vn2r12p"
         ,"Vn2r13p","Vn2r14p","Vn2r15p","Vn2r16p","Vn2r17p","Vn2r18p","Vn2r19p","Vn2r1p","Vn2r20p","Vn2r21p","Vn2r2p","Vn2r3p","Vn2r6p","Vn2r7p","Vn2r9p","Xcr1")

Neuropeptides<-c("Penk","Pomc","Pdyn","Pnoc","Avp","Oxt","Gast","Cck","Sst","Cort","Npvf","Npff","Npy","Ppy","Pyy","Prlh","Calca","Calcb","Iapp","Adm","Adm2","Nppa","Nppb","Nppc","Grp","Nmb","Edn1","Edn2","Edn3","Cgc","Sct","Vip","Adcyap1","Ghrh","Gip","Crh","Ucn","Ucn2","Ucn3","Uts2","Uts2b","Tac1","Tac2","Nms","Nmu","Kng1","Agt","Nts","Chga","Chgb","Scg2","Scg3","Scg5","Vgf","Mln","Ghrl","Gal","Galp","Gnrh1","Npb","Npw","Nps","Nxph1","Nxph2","Nxph3","Nxph4","Ins","Ins1","Ins2","Igf1","Igf2","Rln1","Rln3","Trh","Pthlh","Pmch","Hcrt","Cartpt","Agrp","Prl","Apln","Kiss1","Dbi","Cbln1","Cbln2","Cbln3","Cbln4","Lep","Adipoq","Nampt","Retn","Retnla","Retnlb","Retnlg","Nucb2","Ubl5")

Na_channels_Mouse<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a1","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Scnm1","Scn7a")

K_channels_Mouse<-c("Hcn1","Hcn2","Hcn3","Hcn4","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne2","Kcne3","Kcne4","Kcne1l","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4","Dpp6","Dpp10","Kcnrg","Kcne1l")

Voltage_gated_channels_HUGO<-c("Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Cnga1","Cnga2","Cnga3","Cnga4","Cngb1","Cngb3","Hcn1","Hcn2","Hcn3","Hcn4","Hvcn1","Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Trpa1","Trpc1","Trpc2","Trpc3","Trpc4","Trpc5","Trpc6","Trpc7","Mcoln1","Mcoln2","Mcoln3","Trpm1","Trpm2","Trpm3","Trpm4","Trpm5","Trpm6","Trpm7","Trpm8","Pkd2","Pkd2l1","Pkd2l2","Trpv1","Trpv2","Trpv3","Trpv4","Trpv5","Trpv6","Tpcn1","Tpcn2")

Ligand_gated_ion_channels_HUGO<-c("Htr3a","Htr3b","Htr3c","Htr3d","Htr3e","Cftr","Chrna1","Chrna2","Chrna3","Chrna4","Chrna5","Chrna6","Chrna7","Chrna9","Chrna10","Chrnb1","Chrnb2","Chrnb3","Chrnb4","Chrnd","Chrne","Chrng","Gabra1","Gabra2","Gabra3","Gabra4","Gabra5","Gabra6","Gabrb1","Gabrb2","Gabrb3","Gabrd","Gabre","Gabrg1","Gabrg2","Gabrg3","Gabrp","Gabrq","Gabrr1","Gabrr2","Gabrr3","Gria1","Gria2","Gria3","Gria4","Grid1","Grid2","Grik1","Grik2","Grik3","Grik4","Grik5","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grin3a","Grin3b","Glra1","Glra2","Glra3","Glra4","Glrb","Itpr1","Itpr2","Itpr3","P2rx1","P2rx2","P2rx3","P2rx4","P2rx5","P2rx6","P2rx7","Ryr1","Ryr2","Ryr3","Zacn")

Gray.activity.genes<-c("Arc","Acan","Baiap2","Bdnf","Cdkn1a","Rgs2","Nr4a2","Grasp","Nefm","Nrn1","Ptgs2","Scg2","Vgf","Mest","Bcor","Dusp6","Egr3","Fosb","Fosl2","Nptx2","Actn1","Car12","Dnajb5","Fmnl1","Igsf9b","Inhba","Kcnj2","Lingo1","Lingo3","Mapk4","Pim1","Plk2","Pmepa1","Rtn4rl2","Zdbf2","tyro3","Slitrk3","Aloxe3","Dusp14","Fstl4","Peg10","Sorcs3","Npas4","Nr4a1","Fos","Junb","Egr2","Pcsk1","Egr1","Klf4","Gadd45b","Nr4a3","Dusp1","Fbxo33","Gadd45g","Coq10b","Crem","Btg2","Gla","Trib1","Nfil3","Erf","Per1","Ifrd1","Gpr3","Siah2","Gem","Plk3","Arl5b","Homer1","Eprs","Errfi1","Sstr2","Pcdh8","Stk40","Pvr","Per2","Spry2","Arid5a","Baz1a","Dot1l","Sv2c","Nedd9","Trip10","Etv3","Tgfb2","Pak6","Frmd6","Kcna1","Gcnt2","Etv5","Lonrf3")

KChannels<-c("Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2")
KChannel_reg_subunits<-c("Dpp6","Dpp10","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne1b","Kcne2","Kcne3","Kcne4","Kcne5","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4", "Hcn1","Hcn2","Hcn3","Hcn4")
NaChannels<-c("Asic1","Asic2","Asic3","Asic4","Asic5","Scnn1a","Scnn1b","Scnn1d","Scnn1g","Nalcn","Scn1a","Scn2a","Scn3a","Scn4a","Scn5a","Scn8a","Scn9a","Scn10a","Scn11a","Scn1b","Scn2b","Scn3b","Scn4b","Scnm1")
CaChannels<-c("Cacna1a","Cacna1b","Cacna1c","Cacna1d","Cacna1e","Cacna1f","Cacna1g","Cacna1h","Cacna1i","Cacna1s","Cacna2d1","Cacna2d2","Cacna2d3","Cacna2d4","Cacnb1","Cacnb2","Cacnb3","Cacnb4","Cacng1","Cacng2","Cacng3","Cacng4","Cacng5","Cacng6","Cacng7","Cacng8","Catsperb","Catsperd","Catsperg","Catsper1","Catsper2","Catsper3","Catsper4","Itpr1","Itpr2","Itpr3","Ryr1","Ryr2","Ryr3","Tpcn1","Tpcn2")
ClChannels<-c("Ano1","Ano2","Ano3","Ano4","Ano5","Ano6","Ano7","Ano8","Ano9","Ano10","Best1","Best2","Best3","Best4","Cftr","Clic1","Clic2","Clic3","Clic4","Clic5","Clic6","Clcnka","Clcnkb","Bsnd","Clcn1","Clcn2","Clcn3","Clcn4","Clcn5","Clcn6","Clcn7")

Ion_channels<-c(KChannels,KChannel_reg_subunits,NaChannels,CaChannels,ClChannels)

LTP_genes<-c("Adcy1","Adcy8","Araf","Atf4","Braf","Cacna1c","Calm1","Calm2","Calm3","Calml3","Calml5","Calml6","Camk2a","Camk2b","Camk2d","Camk2g","Camk4","Chp","Chp2","Crebbp","Ep300","Gnaq","Gria1","Gria2","Grin1","Grin2a","Grin2b","Grin2c","Grin2d","Grm1","Grm5","Hras","Itpr1","Itpr2","Itpr3","Kras","Map2k1","Map2k2","Mapk1","Mapk3","Nras","Plcb1","Plcb2","Plcb3","Plcb4","Ppp1ca","Ppp1cb","Ppp1cc","Ppp1r12a","Ppp1r1a","Ppp3ca","Ppp3cb","Ppp3cc","Ppp3r1","Ppp3r2","Prkaca","Prkacb","Prkacg","Prkca","Prkcb","Prkcg","Prkx","Raf1","Rap1a","Rap1b","Rapgef3","Rps6ka1","Rps6ka2","Rps6ka3","Rps6ka6","Braf","Crh","Crhr2","Drd1","Gfap","Gip","Grin2a","Itpr3","Lrrtm1","Lrrtm2","Mapk1","Mecp2","Nlgn1","Nlgn3","Nptn","Nr2e1","Ntrk2","Plk2","Prkcz","Pten","Ptk2b","Ptn","Reln","Rgs14","Rims1","S100b","Serpine2","Shank2","Slc24a2","Slc8a2","Slc8a3","Snap25","Snap47","Snca","Stx3","Stx4","Syt12","Tnr","Vamp2","Creb1","Crtc1","Cx3cr1","Drd1","Drd2","Eif2ak4","Grin2a","Hnrnpk","Nf1","Nlgn3","Nos1","Nptn","Nrgn","Paip2","Ppp1r9a","Ptpn5","Reln","Shank3","Stau1")
Axon_genes<-c("Aak1","Acadm","Adam21","Adam22","Adc","Adcy10","Adcy9","Adcyap1","Adnp","Adora1","Adra2c","Adrb2","Ager","Alcam","Als2","Amigo1","Ank1","Ank3","Anxa3","Anxa5","Ap1s1","Ap2m1","Ap3b1","Ap3b2","Ap3d1","Ap3m1","Ap3m2","Ap3s1","Ap3s2","Apbb1","Apc","App","Atat1","Atcay","Atl1","Atp1a3","Atp6v0d1","Aurka","Avil","Bace1","Bin1","Bloc1s1","Bloc1s2","Bloc1s3","Boc","Bsn","C19orf20","C1orf130","C1orf96","C4orf49","Ca2","Cabp4","Cad","Cadm2","Calb1","Calb2","Calca","Camk2d","Canx","Cck","Cckar","Ccl2","Cdh8","Cdk5","Cdk5r1","Chrm1","Chrm2","Chrm3","Chrm4","Chrna10","Cib1","Clasp2","Cno","Cnr1","Cntf","Cntn2","Cntn4","Cntnap1","Cntnap2","Cobl","Comt","Coro1a","Cplx1","Cplx2","Cpne6","Cpt1c","Creb1","Crh","Crhbp","Crhr2","Ctnna2","Cyp17a1","Dab2ip","Dag1","Dagla","Dcc","Ddc","Dgki","Dicer1","Disc1","Dlg1","Dlg2","Dlg4","Dnajc5","Dnm1","Dnm3","Dock7","Dpysl2","Drd1","Drd2","Drd4","Dscam","Dst","Dtna","Dtnbp1","Dvl1","Eea1","Elk1","Epb41l3","Epha4","Epha5","Ephb1","Ephb2","Ermn","Esr1","Fam168b","Fez1","Fgf13","Fkbp15","Fkbp4","Flrt3","Fmr1","Fxr1","Fzd3","Gabbr1","Gabra2","Gabrg2","Gad2","Gap43","Gars","Gc","Gdpd5","Ghrh","Ghrl","Glul","Gnrh1","Got1","Gper","Gpm6a","Gria1","Gria3","Gria4","Grik2","Grik3","Grik5","Grin1","Grm2","Grm3","Grm7","Gsk3b","Gsto1","Hap1","Hcfc1","Hcn1","Hcn3","Hdac6","Hepacam","Hif1a","Hmbs","Hnrnpk","Hnrnpr","Homer1","Hpca","Htr2a","Htr3a","Htt","Igf2bp1","Ighmbp2","Igsf9","Il1r1","Il1rapl1","Ilk","Iqgap1","Irx3","Itga2","Katna1","Katnb1","Kcna1","Kcna2","Kcna3","Kcna4","Kcna6","Kcnab1","Kcnab2","Kcnb1","Kcnc2","Kcnc4","Kcnh1","Kcnip3","Kcnj11","Kcnk2","Kcnq2","Kcnq3","Kiaa1598","Kif13b","Kif1a","Kif1b","Kif20b","Kif3b","Kif4a","Kirrel3","Klhl20","Klhl24","L1cam","Ldlrap1","Limk1","Llgl1","Lmtk3","Lpar3","Lphn1","Lphn3","Lrfn3","Lrp8","Lrrc4b","Lrrc7","Lrrk2","Lrrtm1","Lsm1","Maf1","Mag","Map3k12","Mapk1","Mapk8ip1","Mapk8ip3","Mapt","Marcks","Mbp","Mink1","Mme","Mt3","Mtmr2","Mtpn","Mul1","Muted","Myh10","Myh14","Myo1d","Myoc","Napa","Ncam2","Ncdn","Ncs1","Ndel1","Nefh","Nefl","Nefm","Nek3","Nf1","Nfasc","Nfib","Ngdn","Nmu","Nov","Npff","Npy1r","Nrcam","Nrg1","Nrgn","Nrn1l","Nrp1","Nrtn","Nrxn1","Ntng2","Ntrk1","Ntrk2","Nts","Ntsr1","Odz3","Olfm1","Omp","Opa1","Ophn1","Oprd1","Oprk1","Oxt","P2rx3","P2rx4","Pacsin1","Pafah1b1","Pak1","Palld","Palm","Pard3","Pard6a","Park7","Pawr","Pcp4","Pdyn","Penk","Pfn2","Pink1","Pldn","Pnoc","Polg","Ppt1","Prkcz","Prss12","Psen1","Ptch1","Ptk2b","Ptpn5","Ptprk","Ptprn","Ptprn2","Ptpro","Ptprz1","Pvalb","Pvrl1","Rab11a","Rab21","Rab3a","Rab5a","Rab7a","Rangap1","Rap1gap","Rapgef3","Ret","Rnf40","Rnf6","Robo1","Robo2","Robo3","Rtn4r","Rufy3","Rxra","Sacs","Sarm1","Scn11a","Scn1a","Scn1b","Scn2a","Scn3a","Scn8a","Sdc3","Sema3a","Sema6a","Sept11","Sept5","Sept6","Serpinf1","Setx","Sirt2","Slc17a7","Slc17a8","Slc18a2","Slc1a2","Slc32a1","Slc38a2","Slc38a7","Slc5a7","Slc6a1","Slc6a3","Slc9a6","Smurf1","Snapin","Snca","Sncb","Sncg","Sod1","Spast","Spg11","Spg7","Sphk1","Spock1","Sptbn1","Sptbn4","Srcin1","Sri","Stat1","Stmn2","Stmn3","Stmn4","Stx6","Stxbp1","Sv2a","Syn1","Syngr1","Synj2","Syp","Syt1","Syt7","Tac1","Tanc1","Tbc1d24","Tgfb1","Tgfb2","Th","Tiam1","Tmem57","Tnfrsf1b","Tnfrsf21","Tpx2","Trpv2","Tubb3","Tubb4a","Tulp1","Uchl1","Ucn","Uhmk1","Unc13a","Unc13b","Vamp1","Vamp2","Vps16","Wdr81","Ywhae","Zfyve27","Znf259")

rollShannon<-function(x){
  p<-x
  cummeRbund:::shannon.entropy(p)
}

my_pseudotime_plot_celltype<-function(cds_subset,marker="Cdh13"){
  markerId<-lookupGeneId(dat.filtered,marker)
  tmp<-melt(data.frame("Pseudotime"=seq(0, max(pData(cds_subset)$Pseudotime),length.out = 100),
                     "celltype1"=as.vector(genSmoothCurves(dat.filtered[markerId,pData(dat.filtered)$celltype==1],cores=1, trend_formula="~sm.ns(Pseudotime,df=3)",relative_expr=T,new_data=newdata)),
                     "celltype2"=as.vector(genSmoothCurves(dat.filtered[markerId,pData(dat.filtered)$celltype==2],cores=1, trend_formula="~sm.ns(Pseudotime,df=3)",relative_expr=T,new_data=newdata))),
          id.vars = "Pseudotime")

  p<-ggplot(meltCDS(dat.filtered,geneset=marker)) + 
    geom_jitter(aes(x=Pseudotime,y=value,color=celltype),width=1) + 
    geom_line(aes(x=Pseudotime,y=value,color=variable,group=variable),data=tmp) + ggtitle(marker) +
    scale_color_manual(values=c(celltype_colors,celltype_colors)) + monocle:::monocle_theme_opts()
  p
}