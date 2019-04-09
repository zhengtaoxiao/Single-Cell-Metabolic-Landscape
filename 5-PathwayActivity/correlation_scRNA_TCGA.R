library(scater)
library(stringr)
library(pheatmap)
source("../utils.R")
library(ggrepel)
library(reshape2)
options(stringsAsFactors = FALSE)


##read the TCGA data
gmt_file <- "../Data/KEGG_metabolism.gmt"

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

#tumors
metabolic_tcga_tumor_sce <- metabolic_tcga_sce[,metabolic_tcga_sce$Sample.Type=="tumor"]
metabolic_tcga_tumor_fpkm <- assay(metabolic_tcga_tumor_sce,"exprs")
################################################
#read the scRNA data
tumor <- "head_neck"
selected_impute_sce <- readRDS(file.path("../2-Imputation/dataset",tumor,"selected_impute_sce.rds"))
normalization_method <- "Deconvolution"
norm_rds_file <- file.path("../3-Normalization/dataset/",tumor,paste0(normalization_method,"_tpm.rds"))
norm_tpm <- readRDS(norm_rds_file)

all_celltypes <- selected_impute_sce$cellType
all_tumors <- as.vector(selected_impute_sce$tumor)
tumors <- unique(all_tumors)
#statistcs of cell numbers for each tumor
for(i in tumors){
  print(i)
  print(table(selected_impute_sce[,selected_impute_sce$tumor==i]$cellType))
}

selected_tumors <- c("HNS17","HNS16","HNS15","HNS12","HNS10","HNS11","HNS8","HNS9","HNS2","HNS1")
corr_list <- list()
cor_pvalues <- list()
for(i in selected_tumors){
  norm_tpm_each <- norm_tpm[,all_tumors==i]
  norm_tpm_each_means <- rowMeans(norm_tpm_each,na.rm = T)
  common_genes <- intersect(names(norm_tpm_each_means),rownames(metabolic_tcga_tumor_fpkm))
  norm_tpm_each_means_select <- norm_tpm_each_means[common_genes]
  metabolic_tcga_tumor_fpkm_select <- metabolic_tcga_tumor_fpkm[common_genes,]
  corrs <- apply(metabolic_tcga_tumor_fpkm_select,2,function(x) cor(x,norm_tpm_each_means_select,use="complete.obs",method="spearman"))
  corrs_pvs <- apply(metabolic_tcga_tumor_fpkm_select,2,function(x) cor.test(x,norm_tpm_each_means_select,use="complete.obs",method="spearman")$p.value)
  corr_list[[i]] <- corrs
  cor_pvalues[[i]] <- -log10(corrs_pvs)
  #write(corrs,"corr.txt",sep=",",ncolumns=length(corrs),append=T)
}


##randomly selected some Tcells to reconstruct bulk
set.seed(1234)
norm_tpm_random <- norm_tpm[,all_celltypes=="T cell"]
norm_tpm_random <- norm_tpm_random[,sample(1:ncol(norm_tpm_random),500)]
norm_tpm_random_mean <- rowMeans(norm_tpm_random)
common_genes <- intersect(names(norm_tpm_random_mean),rownames(metabolic_tcga_tumor_fpkm))
norm_tpm_random_mean_select <- norm_tpm_random_mean[common_genes]
metabolic_tcga_tumor_fpkm_select <- metabolic_tcga_tumor_fpkm[common_genes,]
corrs <- apply(metabolic_tcga_tumor_fpkm_select,2,function(x) cor(x,norm_tpm_random_mean_select,use="complete.obs",method="spearman"))
corrs_pvs <- apply(metabolic_tcga_tumor_fpkm_select,2,function(x) cor.test(x,norm_tpm_random_mean_select,use="complete.obs",method="spearman")$p.value)

