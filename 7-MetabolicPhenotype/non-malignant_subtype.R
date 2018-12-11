library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(pheatmap)
source("../utils.R")
source("runGSEA.R")

pathway_file <- "../Data/KEGG_metabolism.gmt"
hall_gmt <- '../Data/h.all.v6.1.symbols.gmt'

args <- commandArgs()
tumor <- args[6]
outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir) ) dir.create(outDir,recursive=TRUE)

#1. Loading the data
selected_sce <- readRDS(file.path("../1-ReadData/dataset/",tumor,"selected_sce.rds"))
selected_nontumor_sce <- selected_sce[,selected_sce$cellType!="Malignant"]
selected_nontumor_metabolic_sce <- selected_nontumor_sce[rowData(selected_nontumor_sce)$metabolic,]
#=========================================================================
celltypes <- unique(selected_nontumor_metabolic_sce$cellType)

#2. T cell type
tcell_sce <- selected_nontumor_sce[,selected_nontumor_sce$cellType=="T cell"]
tcell_exp <- assay(tcell_sce,"exprs")
tcell_metabolic_sce <- selected_nontumor_metabolic_sce[,selected_nontumor_metabolic_sce$cellType=="T cell"]
tcell_metabolic_exp <- assay(tcell_metabolic_sce,"exprs")

#difference between CD4+ and CD8+
#scatter CD4 and CD8
pdf(file.path(outDir,"CD4_CD8_scatterplot.pdf"))
plot(t(tcell_exp[c("CD4","CD8A"),]))
dev.off()
#select the cells
col_annot <- data.frame(type=rep(NA,ncol(tcell_sce)),row.names=colnames(tcell_sce))
col_annot[(tcell_exp["CD4",]>1)&(tcell_exp["CD8A",]<1),] <- "CD4"
col_annot[(tcell_exp["CD8A",]>1)&(tcell_exp["CD4",]<1),] <- "CD8"

select_tcell_sce <- tcell_sce[,!is.na(col_annot$type)]
select_tcell_metabolic_sce <- tcell_metabolic_sce[,!is.na(col_annot$type)]
select_col_annot <- col_annot[!is.na(col_annot$type),,drop=F]
select_tcell_metabolic_sce$type <- select_col_annot$type
select_tcell_sce$type <- select_col_annot$type
runGSEA(select_tcell_metabolic_sce,"type","CD4","CD8","t",pathway_file,file.path(outDir,"CD4_CD8_GSEA"))

#select CD4 cells, and sort to Tregs and Th cells
cd4_tcell_sce <- select_tcell_sce[,select_tcell_sce$type=="CD4"]
cd4_tcell_exp <- assay(cd4_tcell_sce,"exprs")
tmp <- cd4_tcell_exp[c("FOXP3","IL2RA"),]
cd4_col_annot <- data.frame(type=rep(NA,ncol(cd4_tcell_exp)),row.names=colnames(cd4_tcell_exp))
cd4_col_annot[colSums(tmp)>=2,] <- "Tregs"
cd4_col_annot[colSums(tmp)==0,] <- "Th"
select_cd4_tcell_sce <- cd4_tcell_sce[,!is.na(cd4_col_annot$type)]
select_cd4_tcell_metabolic_sce <- select_cd4_tcell_sce[rowData(select_cd4_tcell_sce)$metabolic,]
cd4_col_annot <- cd4_col_annot[!is.na(cd4_col_annot$type),,drop=F]
select_cd4_tcell_metabolic_sce$type <- cd4_col_annot$type
#plot selection
tmp <- cd4_tcell_exp[c("FOXP3","IL2RA"),]
pdf(file.path(outDir,"FOXP3_CD25_heatmap.pdf"),width=3,height=1)
pheatmap(tmp[,order(colSums(tmp))],show_rownames =T,show_colnames = F,cluster_rows = F,cluster_cols = F)
dev.off()
runGSEA(select_cd4_tcell_metabolic_sce,"type","Tregs","Th","t",pathway_file,file.path(outDir,"Tregs_Ths_GSEA"))


#3 Fibroblast cells: only for head and neck tumors
if (tumor == "head_neck"){
fib_sce <- selected_nontumor_sce[,selected_nontumor_sce$cellType=="Fibroblast"]
fib_exp <- assay(fib_sce,"exprs")
#filter
select <- (fib_exp["FOS", ] >= 1) & (fib_exp["VIM",] >= 1)
select_fib_sce <- fib_sce[,select]
select_fib_exp <- assay(select_fib_sce,"exprs")
#select CAF and Myofib, at least 2 of markers > 1
myofib_markers <- c("ACTA2","MCAM","MYLK","MYL9","PDGFA")
CAFs_markers <- c("FAP","THY1","PDPN","PDGFRA","PDGFRL","MMP2")
select <- (apply(select_fib_exp[myofib_markers,],2,function(x) sum(x >= 1)>=2)) | (apply(select_fib_exp[CAFs_markers,],2,function(x) sum(x >= 1)>=2))
select_fib_sce2 <- select_fib_sce[,select]
select_fib_exp2 <- assay(select_fib_sce2,"exprs")
#write the marker gene
dat <- select_fib_exp2[c("FOS","VIM",myofib_markers,CAFs_markers),]
hr <- hclust(as.dist(1-cor(t(dat),method="pearson")),method="ward.D2")
hc <- hclust(as.dist(1-cor(dat,method="pearson")),method="ward.D2")
mybreaks <- c(seq(-2,0,length.out=ceiling(200/2)+1),
              seq(2/200,2,length.out = floor(200/2)))
library(RColorBrewer)
mycolor=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(200)
pdf(file.path(outDir,"fibro_markers.pdf"),width=3.5,height=2,onefile = T)
pheatmap(dat,show_colnames = F,scale="row",cluster_cols = hc,cluster_rows = hr,breaks=mybreaks,legend = F,color=mycolor)
dev.off()

#cluster to CAFs or myofibroblasts
tmp=select_fib_exp2[c(CAFs_markers,myofib_markers),]
kmeans_res <- kmeans(t(tmp),centers=2)
col_annot <- data.frame(type=rep(NA,ncol(select_fib_exp2)),row.names=colnames(select_fib_exp2))
if(sum(tmp[CAFs_markers,kmeans_res$cluster==1]) > sum(tmp[CAFs_markers,kmeans_res$cluster==2])){
  col_annot[kmeans_res$cluster==1,] <- "CAF"
  col_annot[kmeans_res$cluster==2,] <- "Myofib"
}else{
  col_annot[kmeans_res$cluster==1,] <- "Myofib"
  col_annot[kmeans_res$cluster==2,] <- "CAF"
}
select_fib_metabolic_sce2 <- select_fib_sce2[rowData(select_fib_sce2)$metabolic,]
select_fib_metabolic_sce2$type <- col_annot$type

runGSEA(select_fib_metabolic_sce2,"type","CAF","Myofib","t",pathway_file,file.path(outDir,"CAF_Myofib_GSEA"))
}
date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)