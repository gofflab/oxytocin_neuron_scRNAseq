library(monocle)
library(reshape2)
library(tidyverse)
library(pheatmap)
library(marray)
library(ggrepel)
library(ggbiplot)
library(gridExtra)
library(viridis)
library(VGAM)
library(gplots)


#############
#slackr
#############
library(slackr)

slackrSetup(channel="#dolen_oxytocin_sc", 
            incoming_webhook_url="https://hooks.slack.com/services/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",
            api_token="xoxp-9317854839-9317673428-431659063763-79b68d407c25c45432b5cdb2c8a9ee68") #This is my personal token
############
# Color Palette
############
celltype_colors<-c("darkgreen","gold")
label_colors<-c("darkgreen","gold")

gg_color_hue <- function(n) {
  hues = seq(15, 375, length=n+1)
  hcl(h=hues, l=65, c=100)[1:n]
}

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

myHeatmap<-function(cds,geneset,logMode=FALSE){
  sub<-cds[lookupGeneId(cds,geneset),]
  if(logMode){
    heatmap.2(log10(as.matrix(exprs(sub))+1),scale="row",trace="none",labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(2)[pData(sub)$CellType],labCol=FALSE,distfun=dist,hclustfun=hclust2)
  }else{
    heatmap.2(as.matrix(exprs(sub)),scale="row",trace="none",labRow=fData(sub)$gene_short_name,ColSideColors=gg_color_hue(2)[pData(sub)$CellType],labCol=FALSE,distfun=dist,hclustfun=hclust2)
  }
}

myPHeatmap<-function(cds,geneset,logMode=TRUE){
  sub<-cds[lookupGeneId(cds,geneset),]
  mat<-as.matrix(exprs(sub))
  if(logMode){
    mat<-log10(mat+1)
  }
  pheatmap(mat=log10(as.matrix(exprs(sub))+1),
            scale="row",
            labels_row=fData(sub)$gene_short_name,
            annotation_col=pData(sub)[,c("Cell_Type_guess","source_plate","Cluster","Fluorogold","Total_mRNAs","num_genes_expressed")],
            labelsCol=FALSE)
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
  dat.1<-apply(dat[gene_ids,pData(cds)$Cluster %in% cluster_1],1,mean)
  dat.2<-apply(dat[gene_ids,pData(cds)$cCluster %in% cluster_2],1,mean)
  total.diff<-abs(dat.1-dat.2)
  genes.diff<-names(which(total.diff>thresh.use))
  genes.use<-genes.diff[genes.diff%in%rownames(dat)]
  res<-per_gene_AUC(dat[,pData(cds)$Cluster %in% cluster_1],dat[,pData(cds)$Cluster %in% cluster_2],genes.use)
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

#############
# Relevant Gene Sets
#############

KChannels<-c("Kcnma1","Kcnn1","Kcnn2","Kcnn3","Kcnn4","Kcnu1","Kcnt1","Kcnt2","Kcnk1","Kcnk2","Kcnk3","Kcnk4","Kcnk5","Kcnk6","Kcnk7","Kcnk9","Kcnk10","Kcnk12","Kcnk13","Kcnk15","Kcnk16","Kcnk17","Kcnk18","Kcnj1","Kcnj2","Kcnj3","Kcnj4","Kcnj5","Kcnj6","Kcnj8","Kcnj9","Kcnj10","Kcnj11","Kcnj12","Kcnj13","Kcnj14","Kcnj15","Kcnj16","Kcnj18","Kcna1","Kcna2","Kcna3","Kcna4","Kcna5","Kcna6","Kcna7","Kcna10","Kcnb1","Kcnb2","Kcnc1","Kcnc2","Kcnc3","Kcnc4","Kcnd1","Kcnd2","Kcnd3","Kcnf1","Kcng1","Kcng2","Kcng3","Kcng4","Kcnh1","Kcnh2","Kcnh3","Kcnh4","Kcnh5","Kcnh6","Kcnh7","Kcnh8","Kcnq1","Kcnq2","Kcnq3","Kcnq4","Kcnq5","Kcns1","Kcns2","Kcns3","Kcnv1","Kcnv2")
KChannel_reg_subunits<-c("Dpp6","Dpp10","Kcnab1","Kcnab2","Kcnab3","Kcne1","Kcne1b","Kcne2","Kcne3","Kcne4","Kcne5","Kcnip1","Kcnip2","Kcnip3","Kcnip4","Kcnmb1","Kcnmb2","Kcnmb3","Kcnmb4", "Hcn1","Hcn2","Hcn3","Hcn4")
