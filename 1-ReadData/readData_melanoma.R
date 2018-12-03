source("../utils.R")
library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(reshape2)
outdir <- "dataset/melanoma"
if(!dir.exists(outdir)) dir.create(outdir)
#################################################################################################
# 1. Read the data
#################################################################################################
raw_tpm_file <- "dataset/GSE72056_melanoma_single_cell_corrected.txt"
tmp_data <- read.table(raw_tpm_file,head=T,sep="\t",quote=NULL,stringsAsFactors=F)

#there are two duplicated genes: MARCH1 and MARCH2
#For each duplicated gene, remove the one with lower mean expression
remove_rows <- c()
march1_rows <- grep("^MARCH1$",tmp_data$Cell)
march1_rows_mean <- apply(tmp_data[march1_rows,2:ncol(tmp_data)],1,mean)
remove_rows <- c(remove_rows,names(march1_rows_mean)[order(march1_rows_mean)][1])

march2_rows <- grep("^MARCH2$",tmp_data$Cell)
march2_rows_mean <- apply(tmp_data[march2_rows,2:ncol(tmp_data)],1,mean)
remove_rows <- c(remove_rows,names(march2_rows_mean)[order(march2_rows_mean)][1])
remove_rows <- as.integer(remove_rows)


#column data
col_data = data.frame(tumor=t(tmp_data[1,2:ncol(tmp_data)]),malignant=t(tmp_data[2,2:ncol(tmp_data)]),
                      cellType=t(tmp_data[3,2:ncol(tmp_data)]))
colnames(col_data) <- c("tumor","malignant","cellType")
col_data$tumor <- factor(paste0("T",col_data$tumor))
#malignant(1=no,2=yes,0=unresolved)
#non-malignant cell type (0=unclassfied,1=T,2=B,3=Macro.4=Endo.,5=CAF;6=NK)
col_data[col_data$malignant==1,"malignant"] <- "Non-malignant"
col_data[col_data$malignant==2,"malignant"] <- "Malignant"
col_data[col_data$malignant==0,"malignant"] <- "Unresolved"
col_data$malignant <- factor(col_data$malignant,levels=c("Non-malignant","Malignant","Unresolved"))
col_data[col_data$cellType==0,"cellType"] <- "Unknow"
col_data[col_data$cellType==1,"cellType"] <- "T cell"
col_data[col_data$cellType==2,"cellType"] <- "B cell"
col_data[col_data$cellType==3,"cellType"] <- "Macrophage"
col_data[col_data$cellType==4,"cellType"] <- "Endothelial"
col_data[col_data$cellType==5,"cellType"] <- "CAF"
col_data[col_data$cellType==6,"cellType"] <- "NK"
#some unknow is malignant
tumor_select <- (col_data$cellType=="Unknow") & (col_data$malignant=="Malignant")
col_data[tumor_select,"cellType"] <- "Malignant"
col_data$cellType <- factor(col_data$cellType,levels = c("Unknow","Malignant","T cell","B cell","Macrophage","Endothelial","CAF","NK"))

#remove the annotation lines
remove_rows <- c(remove_rows,1,2,3)
all_data <- tmp_data[-remove_rows,]
rownames(all_data) <- all_data$Cell
all_data$Cell=NULL
rm(tmp_data)

####################################################################################################
###2. marker the metabolic genes
####################################################################################################

pathways <- gmtPathways("../Data/KEGG_metabolism.gmt")
metabolics <- unique(as.vector(unname(unlist(pathways))))
row_data <- data.frame(metabolic=rep(FALSE,nrow(all_data)),row.names = rownames(all_data))
row_data[rownames(row_data)%in%metabolics,"metabolic"]=TRUE

####################################################################################################
#3. build scater object
####################################################################################################
all_data <- data.matrix(all_data)
raw_tpm <- 2^ all_data - 1
sce <- SingleCellExperiment(
  assays = list(tpm=raw_tpm,exprs=all_data),
  colData = col_data,
  rowData = row_data
)

####################################################################################################
#4. filtering the cells.
# cutoff 50 
####################################################################################################

#malignant cells
tumor_sce <- sce[,sce$cellType == "Malignant"]
nontumor_sce <- sce[,!sce$cellType%in%c("Unknow","Malignant")]
#select tumor cells
tumor_sample_stats <- table(tumor_sce$tumor)
tumor_sample_select <- names(tumor_sample_stats)[tumor_sample_stats>=50]
selected_tumor_sce <- tumor_sce[,tumor_sce$tumor %in% tumor_sample_select]

#select notumor
nontumor_sample_stats <- table(nontumor_sce$cellType)
nontumor_sample_select <- names(nontumor_sample_stats)[nontumor_sample_stats>=50]
selected_nontumor_sce <- nontumor_sce[,nontumor_sce$cellType %in% nontumor_sample_select]
#selected_sce
selected_columns <- unique(c(colnames(selected_tumor_sce),colnames(selected_nontumor_sce)))
selected_sce <- sce[,colnames(sce) %in% selected_columns]

#rename the patients
unique_tumors <- unique(selected_sce$tumor)
levels(selected_sce$tumor) = paste0("MEL",seq(1,length(unique_tumors)))


###5. Save as R data friles
#saveRDS(sce,file.path(outdir,"sce.rds"))
saveRDS(selected_sce,file.path(outdir,"selected_sce.rds"))
