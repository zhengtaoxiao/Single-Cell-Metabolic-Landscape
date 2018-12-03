source("../utils.R")
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
library(plyr)

outdir <- "dataset/head_neck"
if(!dir.exists(outdir)) dir.create(outdir)
# 1. Read the ata
#################################################################################################
#The value is normalized Expression value
raw_tpm_file <- "dataset/GSE103322_HNSCC_all_data.txt"
tmp_data <- read.table(raw_tpm_file,head=T,sep="\t",row.names=1,quote="\'",stringsAsFactors=F)

tumor <- sapply(str_split(colnames(tmp_data),"_"),function(x) x[1])
tumor <- str_sub(tumor,-2,-1)
tumor <- paste0("MEEI",str_replace(tumor,"C",""))

cell_type <- as.character(tmp_data[5,])
malignant <- as.character(tmp_data[3,]) == "1"   
cell_type[malignant] <- "Malignant"
cell_type[cell_type==0] <- "Unknow" 

#column data
col_data <- data.frame(tumor=tumor,cellType=cell_type,
                       lymph=as.integer(tmp_data[2,]),
                       row.names=colnames(tmp_data))
#remove the annotation lines
remove_rows <- c(1,2,3,4,5)
all_data <- tmp_data[-remove_rows,]
rm(tmp_data)

##############################################################################################
###2. marker the metabolic genes
pathways <- gmtPathways("../Data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE


####################################################################################################
#3. build scater object
all_data <- data.matrix(all_data)
raw_tpm <- 2^ all_data - 1
sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm, exprs=all_data),
  colData = col_data,
  rowData = row_data
)


####################################################################################################
#4. cell types
#malignant cells
tumor_sce <- sce[,sce$cellType == "Malignant"]
nontumor_sce <- sce[,!sce$cellType%in%c("Unknow","Malignant")]
#select tumor cells
tumor_sample_stats <- table(tumor_sce$tumor)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats>=50]
selected_tumor_sce <- tumor_sce[,tumor_sce$tumor %in% tumor_sample_select]
#select notumor
nontumor_stats <- table(nontumor_sce$cellType)
nontumor_select <- names(nontumor_stats)[nontumor_stats>=50]
selected_nontumor_sce <- nontumor_sce[,nontumor_sce$cellType %in% nontumor_select]
#select sce
selected_columns <- unique(c(colnames(selected_tumor_sce),colnames(selected_nontumor_sce)))
selected_sce <- sce[,colnames(sce) %in% selected_columns]

#rename the patients
selected_sce$tumor <- factor(selected_sce$tumor)
selected_sce$tumor <- mapvalues(selected_sce$tumor, from=c("MEEI5","MEEI6","MEEI7","MEEI8","MEEI10","MEEI12","MEEI13","MEEI16","MEEI17","MEEI18","MEEI20","MEEI22","MEEI23","MEEI24","MEEI25","MEEI26","MEEI28","MEEIC"),
                              to=paste0("HNS",seq(18))) 

#saveRDS(sce,file.path(outdir,"sce.rds"))
saveRDS(selected_sce,file.path(outdir,"selected_sce.rds"))
