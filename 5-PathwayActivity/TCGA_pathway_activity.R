library(scater)
library(stringr)
library(pheatmap)
source("../utils.R")
library(ggrepel)
library(reshape2)
library(doParallel)
options(stringsAsFactors = FALSE)

gmt_file <- "../Data/KEGG_metabolism.gmt"
num_cores <- 4

tumor <- "HNSC" 
outDir <- file.path("TCGA",tumor)
if(!dir.exists(outDir)) dir.create(outDir,recursive = T)
TCGA_exp_file <- file.path("../Data/TCGA/",tumor,"exp.txt")
TCGA_type_file <- file.path("../Data/TCGA/",tumor,"type.txt")

tryCatch({
  exp_data <- read.table(TCGA_exp_file,sep="\t",header=TRUE,row.names=1)
  tumor_type <- read.table(TCGA_type_file,sep="\t",header=TRUE,row.names=1)
},error = function(e){
  stop("TCGA expression file can't be opened, please generate it in ../Data/TCGA folder")
})

tumor_type$Sample.Type <- gsub("Primary Tumor","tumor",tumor_type$Sample.Type)
tumor_type$Sample.Type <- gsub("Metastatic","tumor",tumor_type$Sample.Type)
tumor_type$Sample.Type <- gsub("Solid Tissue Normal","normal",tumor_type$Sample.Type)
rownames(tumor_type) <- str_replace_all(rownames(tumor_type),"-",".")
row_data <- data.frame(geneid = str_split_fixed(rownames(exp_data),fixed('.'),2)[,1],
                       row.names=rownames(exp_data))
#get the gene name 
library("biomaRt")
ensembl <- useMart("ensembl",host="http://www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene_names <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                    filters="ensembl_gene_id",
                    values = row_data$geneid,
                    mart = ensembl)
#remove the duplicated"
gene_names <- gene_names[!duplicated(gene_names$ensembl_gene_id),]
rownames(gene_names) <- gene_names$ensembl_gene_id
row_data$genename <- gene_names[row_data$geneid,"hgnc_symbol"]
#build sce
tcga_sce <- SingleCellExperiment(
	assays = list(fpkm=data.matrix(exp_data), exprs=log2(data.matrix(exp_data)+1)),
	colData = tumor_type,
	rowData = row_data
)

### only select the matched samples.
tumor_ids <- rownames(tumor_type)[tumor_type$Sample.Type=="tumor"]
normal_ids <- rownames(tumor_type)[tumor_type$Sample.Type=="normal"]
tumor_ids <- str_sub(tumor_ids,start=1,end=12)
normal_ids <- str_sub(normal_ids,start=1,end=12)
common_ids <- intersect(tumor_ids,normal_ids)
all_ids <- str_sub(rownames(tumor_type),start=1,end=12)
selection = all_ids %in% common_ids
matched_tcga_sce <- tcga_sce[,selection]
#marker the metabolic genes
pathways <- gmtPathways(gmt_file)
metabolics <- unique(as.vector(unlist(pathways)))
rowData(matched_tcga_sce)$metabolic <- FALSE
rowData(matched_tcga_sce)[rowData(matched_tcga_sce)$genename %in% metabolics,"metabolic"] <- TRUE

metabolic_tcga_sce <- matched_tcga_sce[rowData(matched_tcga_sce)$metabolic,]
#change the gene id to genename
rownames(metabolic_tcga_sce) <- rowData(metabolic_tcga_sce)$genename
metabolic_tcga_fpkm <- assay(metabolic_tcga_sce,"fpkm")
metabolic_tcga_fpkm <- metabolic_tcga_fpkm[rowSums(metabolic_tcga_fpkm)>0,]

################################################ pathway activity
pathway_number  <- num_of_pathways(gmt_file,rownames(metabolic_tcga_sce)[rowData(metabolic_tcga_sce)$metabolic])
pathways <- gmtPathways(gmt_file)
pathway_names <- names(pathways)
mean_expression_noshuffle <- matrix(NA,nrow=length(pathway_names),ncol=2,dimnames = list(pathway_names,c("tumor","normal")))
mean_expression_shuffle <- matrix(NA,nrow=length(pathway_names),ncol=2,dimnames = list(pathway_names,c("tumor","normal")))

##Calculate the pathway activities
cl <- makeCluster(num_cores)
registerDoParallel(cl)

for(p in pathway_names){
  genes <- pathways[[p]]
  genes_comm <- intersect(genes, rownames(metabolic_tcga_sce))
  if (length(genes_comm) <5) next
  
  pathway_metabolic_tpm <- assay(metabolic_tcga_sce,"fpkm")[genes_comm,]
  mean_exp_eachCellType <- apply(pathway_metabolic_tpm,1,function(x) by(x,metabolic_tcga_sce$Sample.Type,mean))
  ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
  
  pathway_number_weight = 1 / pathway_number[genes_comm,]
  mean_exp_pathway <- apply(ratio_exp_eachCellType,2,function(x) weighted.mean(x,pathway_number_weight/sum(pathway_number_weight)))

  mean_expression_noshuffle[p,] <-  mean_exp_pathway[c("tumor","normal")]
  
  ###shuffle 
  tmp <- foreach(i=1:1000) %dopar%{
    shuffle_cell_types <- sample(metabolic_tcga_sce$Sample.Type)
    mean_exp_eachCellType <- apply(pathway_metabolic_tpm,1,function(x) by(x,shuffle_cell_types,mean))
    ratio_exp_eachCellType <- t(mean_exp_eachCellType) / colMeans(mean_exp_eachCellType)
    mean_exp_pathway <- apply(ratio_exp_eachCellType,2,function(x) weighted.mean(x,pathway_number_weight/sum(pathway_number_weight)))
    mean_exp_pathway[c("tumor","normal")]
  }
  shuffle_results <- matrix(unlist(tmp), ncol=2,byrow=T)
  rownames(shuffle_results) <- 1:1000
  colnames(shuffle_results) <- c("tumor","normal")
  for(c in c("tumor","normal")){
    if(is.na(mean_expression_noshuffle[p,c])) next
    if(mean_expression_noshuffle[p,c]>1){
      pval <- (1+sum(shuffle_results[,c] > mean_expression_noshuffle[p,c])) / 1000
    }else if(mean_expression_noshuffle[p,c]<1){
      pval <- (1+sum(shuffle_results[,c] < mean_expression_noshuffle[p,c])) / 1000 
    }
    mean_expression_shuffle[p,c] <- mean_expression_noshuffle[p,c]
    if(pval>0.01) mean_expression_shuffle[p,c] <- NA
  }
}
stopCluster(cl)
##write
write.table(mean_expression_shuffle,file.path(outDir,"pathway_score_shuffle.txt"),quote=F,sep="\t")
write.table(mean_expression_noshuffle,file.path(outDir,"pathway_score_noshuffle.txt"),quote=F,sep="\t")

not_na <- apply(mean_expression_noshuffle,1,function(x) sum(is.na(x)))
mean_expression_noshuffle <- mean_expression_noshuffle[not_na < 1,]
dat <- mean_expression_noshuffle[rowMaxs(mean_expression_noshuffle)>0,]



#row sort
sort_row <- c()
sort_column <- c()

for(i in colnames(dat)){
  select_row <- which(rowMaxs(dat,na.rm = T) == dat[,i])
  tmp <- rownames(dat)[select_row][order(dat[select_row,i],decreasing = T)]
  sort_row <- c(sort_row,tmp)
}
sort_column <- apply(dat[sort_row,],2,function(x) order(x)[nrow(dat)])
sort_column <- names(sort_column)
dat[is.na(dat)] <- 1

#order according the pathway activity
pdf(file.path(outDir, "KEGGpathway_activity_heatmap.pdf"),onefile=T,width=6,height=9)
mybreaks <- c(
  seq(0, 0.5, length.out=33),
  seq(0.51, 1.5, length.out=33),
  seq(1.51, max(dat),length.out=34)
) 
color <- colorRampPalette(c("blue","white","red"))(100)
pheatmap(dat[sort_row,sort_column],cluster_cols = F,cluster_rows = F,color=color,breaks=mybreaks)
dev.off()

#pathway activity distribution
bulk_df <- melt(mean_expression_noshuffle)
g <- ggplot(bulk_df,aes(x=Var2,y=value,fill=Var2)) +
  geom_violin(trim=F,size=0.2,show.legend = F,width=1.0) + labs(y=NULL,x=NULL) + 
  stat_summary(fun.y = median,geom="point",size=1,color="blue") +
  scale_fill_brewer(palette="Set2") +
  theme_classic() + 
  scale_y_continuous(limits=c(0,3),breaks=0:3,labels=0:3)+
  theme(legend.position="none",
        axis.text.x=element_text(colour="black", size = 6,angle=45,hjust=1,vjust=1),
        axis.text.y=element_text(colour="black", size = 6),
        axis.line=element_line(size=0.2,color="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.ticks.length= unit(.5, "mm"))

ggsave(file.path(outDir,"pathway_activity_violinplot.pdf"),g,width = 0.8,height=1.3,units="in",device="pdf",useDingbats=FALSE)


## correlation iwth scRNA data
sc_dat <- read.table(file.path("dataset/head_neck/KEGGpathway_activity_noshuffle.txt"),header=T,row.names=1,sep="\t")
dat <- cbind(sc_dat[sort_row,"Malignant",drop=T], dat[,"tumor",drop=T])
pdf(file.path(outDir,"bulk_scRNA_pathway_activity_correlation.pdf"))
plot(x=dat[,2],y=dat[,1])
dev.off()

#
