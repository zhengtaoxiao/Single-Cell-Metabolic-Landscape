library(scater)
library(stringr)
library("Rtsne")
library(pheatmap)
library(RColorBrewer)
source("../utils.R")
options(stringsAsFactors=FALSE)


args <- commandArgs()
tumor <- args[6]

outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1. Loading the data
selected_sce <- readRDS(file.path("../2-Imputation/dataset",tumor,"selected_impute_sce.rds"))
selected_tumor_sce <- selected_sce[,selected_sce$cellType == "Malignant"]
selected_tumor_metabolic_sce <- selected_tumor_sce[rowData(selected_tumor_sce)$metabolic,]

#correlation matrix
dat <- assay(selected_tumor_metabolic_sce,"exprs")
dist_dat  <- cor(dat,method="spearman")
hc <- hclust(as.dist(1-dist_dat),method="ward.D2")
mycolor = colorRampPalette(c("white","red"))(10)
mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(8)
col_row <- as.data.frame(colData(selected_tumor_metabolic_sce)[,"tumor",drop=F])
col_row$tumor <- factor(col_row$tumor)
pdf(file.path(outDir,"malignant_metabolic_correlationMatrix.pdf"),width=5,height=4,onefile=T)
pheatmap(dist_dat,show_rownames = F,show_colnames = F,cluster_rows = hc,cluster_cols =hc,annotation_row = col_row,color=mycolor, annotation_legend = T)
dev.off()

####### For melanoma tumor, move MEL3 to after MEL12.
####### 
if(tumor == "melanoma"){
  order_names <- hc$labels[hc$order]
  MEL12_names <- rownames(col_row)[col_row$tumor=="MEL12"]
  order_names2 <- MEL12_names
  MEL3_names <- rownames(col_row)[col_row$tumor=="MEL3"]
  order_names2 <- c(order_names2,MEL3_names)
  other_names <- order_names
  other_names <- order_names[!(order_names %in% MEL3_names)]
  other_names <- other_names[!(other_names %in% MEL12_names)]
  order_names2 <- c(order_names2,other_names)
  dist_dat2 <- dist_dat[order_names2,order_names2]
  pdf(file.path(outDir,"malignant_metabolic_correlationMatrix2.pdf"),width=5,height=4,onefile=T)
  pheatmap(dist_dat2,show_rownames = F,show_colnames = F,cluster_rows = F,cluster_cols =F,annotation_row = col_row,color=mycolor, annotation_legend = T)
  dev.off()
  
}


