library("biomaRt")
library(stringr)
source("../utils.R")

if (!file.exists("CCLE_RNAseq_rsem_genes_tpm_20180929.txt")){
    print("Please download CCLE data from https://portals.broadinstitute.org/ccle or https://figshare.com/articles/scRNA-seq_Datasets/7174922/2")
    q()
}

headLine <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",nrow=1,header=F,sep="\t")
colClass <- rep("numeric",length(headLine))
colClass[1] <- "character"
colClass[2] <- "NULL"
ccld_data <- read.table("CCLE_RNAseq_rsem_genes_tpm_20180929.txt",header=T,sep="\t",row.names=1,colClasses = colClass)
rownames(ccld_data) <- str_split_fixed(rownames(ccld_data),fixed('.'),2)[,1]
ensembl <- useMart("ensembl",host="http://www.ensembl.org")
ensembl <- useDataset("hsapiens_gene_ensembl",mart=ensembl)
gene_names <- getBM(attributes=c("ensembl_gene_id","hgnc_symbol"),
                    filter="ensembl_gene_id",
                    values=rownames(ccld_data),
                    mart=ensembl)
gene_names <- gene_names[!is.null(gene_names[,2]),]
gene_names <- gene_names[!duplicated(gene_names[,2]),]
ccld_data <- ccld_data[gene_names[,1],]
rownames(ccld_data) <- gene_names[,2]

pathway_file <- "../Data/KEGG_metabolism.gmt"
hallmark_gmt <- "../Data/h.all.v6.1.symbols.gmt"

metabolic_pathways <- gmtPathways(pathway_file)
hallmarks <- gmtPathways(hallmark_gmt)

oxphos_genes <- metabolic_pathways[["Oxidative phosphorylation"]]
glycolysis_genes <- metabolic_pathways[["Glycolysis / Gluconeogenesis"]]
hypoxia_genes <- hallmarks[["HALLMARK_HYPOXIA"]]


oxphos_exp <- ccld_data[rownames(ccld_data)%in% oxphos_genes,]
glycolysis_exp <- ccld_data[rownames(ccld_data)%in%glycolysis_genes,]
hypoxia_exp <- ccld_data[rownames(ccld_data)%in%hypoxia_genes,]

oxphos_mean <- colMeans(oxphos_exp,na.rm=T)
glycolysis_mean <- colMeans(glycolysis_exp,na.rm=T)
hypoxia_mean <- colMeans(hypoxia_exp,na.rm=T)

dat <- data.frame(OXPHOS=oxphos_mean,Glycolysis=glycolysis_mean,Hypoxia=hypoxia_mean)
print("correlation:")
print(cor(da))
