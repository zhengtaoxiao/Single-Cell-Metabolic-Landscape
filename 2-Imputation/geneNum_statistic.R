#comparing the metabolic gene number before and after imputation

library(scater)
options(stringsAsFactors=FALSE)

args <- commandArgs()
tumor <- args[6]

outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1 load the data
selected_impute_sce <- readRDS(file.path("./dataset",tumor,"selected_impute_sce.rds"))
selected_sce <- readRDS(file.path("../1-ReadData/dataset/",tumor,"selected_sce.rds"))
cell_types <- unique(selected_sce$cellType)

## check the expressed gene number before and after imputation
gene_number_before <- list()
gene_number_after <- list()

for(c in cell_types){
  each_exp <- assay(selected_sce[rowData(selected_sce)$metabolic,selected_sce$cellType==c],"exprs")
  expressed_genes <- apply(each_exp,2,function(x) sum(x>0))
  impute_each_exp <- assay(selected_impute_sce[rowData(selected_impute_sce)$metabolic,selected_impute_sce$cellType==c],"exprs")
  impute_expressed_genes <- apply(impute_each_exp,2,function(x) sum(x>0))
  gene_number_before[[c]] <- unname(expressed_genes)
  gene_number_after[[c]] <- unname(impute_expressed_genes)
}
pdf(file.path(outDir,"metabolicGeneNumber_beforeImputation.pdf"))
boxplot(gene_number_before,ylim=c(0,1500))
dev.off()
pdf(file.path(outDir,"metabolicGeneNumber_afterImputation.pdf"))
boxplot(gene_number_after)
dev.off()
