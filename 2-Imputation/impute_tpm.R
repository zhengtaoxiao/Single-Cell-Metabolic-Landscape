library(scImpute)
library(scater)

args <- commandArgs()
tumor <- args[6]
num_cores <- 4 #for windows it should be 1

selected_sce <- readRDS(file.path("../1-ReadData/dataset/",tumor,"selected_sce.rds"))

outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir)) dir.create(outDir,recursive=TRUE)

#impute the tumor and non-tumor seperately
selected_tumor_sce <- selected_sce[,selected_sce$cellType=="Malignant"]
selected_nontumor_sce <- selected_sce[,selected_sce$cellType!="Malignant"]

#write the tpm matrix
selected_tumor_tpm <- tpm(selected_tumor_sce) 
selected_nontumor_tpm <- tpm(selected_nontumor_sce) 
labels_tumor <- selected_tumor_sce$tumor
labels_nontumor <- selected_nontumor_sce$cellType

write.csv(selected_tumor_tpm,file.path(outDir,"tumor.tpm"))
write.csv(selected_nontumor_tpm,file.path(outDir,"nontumor.tpm"))

##prepare the gene length file
all_gene_lengths <- read.table("../Data/gene_length.txt",sep="\t",header=F,row.names=1)
tmp <- intersect(rownames(all_gene_lengths),rownames(selected_tumor_tpm))
if (length(tmp) != nrow(selected_tumor_tpm)){
  warning("check the length file")
	print(setdiff(rownames(selected_tumor_tpm),rownames(all_gene_lengths)))
	q()
}
genelen <- all_gene_lengths[rownames(selected_tumor_tpm),]
genelen <- as.numeric(as.vector(genelen))
scimpute(file.path(outDir,"tumor.tpm"),infile="csv",outfile="csv",out_dir=file.path(outDir,"malignant_"),
	labeled=TRUE,labels=as.vector(labels_tumor),
	type="TPM",genelen=genelen,drop_thre=0.5,ncores=num_cores)

imputed_tpm <- read.csv(file.path(outDir,"malignant_scimpute_count.csv"),header=T,row.names=1)
tpm(selected_tumor_sce) <- data.matrix(imputed_tpm) 
assay(selected_tumor_sce,"exprs") <- data.matrix(log2(imputed_tpm + 1))

#non tumor
scimpute(file.path(outDir,"nontumor.tpm"),infile="csv",outfile="csv",out_dir=file.path(outDir,"non-malignant_"),
	labeled=TRUE,labels=as.vector(labels_nontumor),
	type="TPM",genelen=genelen,drop_thre=0.5,ncores=num_cores)

imputed_tpm <- read.csv(file.path(outDir,"non-malignant_scimpute_count.csv"),header=T,row.names=1)
tpm(selected_nontumor_sce) <- data.matrix(imputed_tpm) 
assay(selected_nontumor_sce,"exprs") <- data.matrix(log2(imputed_tpm + 1))

#save as sce
impute_tpm <- cbind(tpm(selected_tumor_sce), tpm(selected_nontumor_sce))
impute_exprs <- cbind(assay(selected_tumor_sce,"exprs"),assay(selected_nontumor_sce,"exprs"))
impute_tpm <- impute_tpm[,colnames(selected_sce)]
impute_exprs <- impute_exprs[,colnames(selected_sce)]

selected_impute_sce <- SingleCellExperiment(
  assays = list(tpm = impute_tpm, exprs=impute_exprs),
  colData = colData(selected_sce),
  rowData = rowData(selected_sce)
)
saveRDS(selected_impute_sce,file.path(outDir,"selected_impute_sce.rds"))
