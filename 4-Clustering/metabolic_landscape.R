library(scater)
library(reshape2)
options(stringsAsFactors=FALSE)

args <- commandArgs()
tumor <- args[6]
outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#1. Loading the tumor data
selected_sce <- readRDS(file.path("../2-Imputation/dataset",tumor,"selected_impute_sce.rds"))

selected_metabolic_sce <- selected_sce[rowData(selected_sce)$metabolic,]
selected_tumor_metabolic_sce <- selected_metabolic_sce[,selected_metabolic_sce$cellType=="Malignant"]
selected_nontumor_metabolic_sce <- selected_metabolic_sce[,selected_metabolic_sce$cellType!="Malignant"]

#===========================================================================
# tsne for tumor cells
set.seed(12345)
library("Rtsne")
tsne_out <- Rtsne(t(assay(selected_tumor_metabolic_sce,"exprs")),initial_dims=20,theta=0.0,perplexity = 30)
tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_tumor_metabolic_sce)$tumor)
g <- ggplot(tmp) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
ggsave(file.path(outDir,"tumor_metabolic_tsne.pdf"),g,width=4,height=3)

#===========================================================================
# tsne for nontumors
tsne_out <- Rtsne(t(assay(selected_nontumor_metabolic_sce,"exprs")),initial_dims=20,perplexity=30,theta=0.0)
tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_nontumor_metabolic_sce)$tumor)

g <- ggplot(tmp) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
ggsave(file.path(outDir,"nontumor_metabolic_tsne_tumorColor.pdf"),g,width=4,height=3)
tmp <- data.frame(x=tsne_out$Y[,1],y=tsne_out$Y[,2],group = colData(selected_nontumor_metabolic_sce)$cellType)
g <- ggplot(tmp) + geom_point(aes(x, y, colour = group), size = 1) +
  labs(x = "tSNE1",y = "tSNE2") +theme_bw() + ggtitle("Rtsne")
ggsave(file.path(outDir,"nontumor_metabolic_tsne_cellTypeColor.pdf"),g,width=4,height=3)
