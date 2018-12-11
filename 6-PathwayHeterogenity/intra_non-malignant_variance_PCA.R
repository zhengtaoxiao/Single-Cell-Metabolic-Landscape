library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(pheatmap)
library(gtools)
library(ggplot2)
source("../utils.R")
source("runGSEA_preRank.R")

args <- commandArgs()
tumor <- args[6]
outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir) ) dir.create(outDir,recursive=TRUE)
pathway_file <- "../Data/KEGG_metabolism.gmt"

#1. Loading the data
selected_sce <- readRDS(file.path("../1-ReadData/dataset/",tumor,"selected_sce.rds"))
selected_nontumor_sce <- selected_sce[,selected_sce$cellType!="Malignant"]
selected_nontumor_metabolic_sce <- selected_nontumor_sce[rowData(selected_nontumor_sce)$metabolic,]
#=========================================================================
celltypes <- unique(selected_nontumor_metabolic_sce$cellType)

#2.Tumor cells
enrich_data_df <- data.frame(x=NULL,y=NULL,NES=NULL,PVAL=NULL)
pc_plotdata <- data.frame(x=numeric(),y=numeric(),
                          sel=character(),types=character())

for (t in celltypes){
  t2 <- str_replace(t," ","")
  each_metabolic_sce <- selected_nontumor_metabolic_sce[,selected_nontumor_metabolic_sce$cellType==t]
  each_metabolic_tpm <- assay(each_metabolic_sce,"exprs")
  each_metabolic_tpm <- each_metabolic_tpm[rowSums(each_metabolic_tpm)>0,]

  x <- each_metabolic_tpm
  ntop <- nrow(x)
  rv <- rowVars(x)
  select <- order(rv, decreasing=TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(x[select,]))
  percentVar <- pca$sdev^2 / sum( pca$sdev^2 )
  
  ###select PCs that explain at least 80% of the variance
  cum_var <- cumsum(percentVar)
  select_pcs <- which(cum_var>0.8)[1]
  ###plot the PCA and explained variances
  tmp_plotdata <- data.frame(x=1:length(percentVar),y=percentVar,
                             sel=c(rep("y",select_pcs),rep("n",length(percentVar)-select_pcs)),
                             types=rep(t,length(percentVar)))
  pc_plotdata <- rbind(pc_plotdata,tmp_plotdata)
  ###
  
  pre_rank_matrix <- as.matrix(rowSums(abs(pca$rotation[,1:select_pcs])))
  runGSEA_preRank(pre_rank_matrix,pathway_file,t2)

  #get the result
  result_dir <- list.files(path="preRankResults",pattern = paste0("^",t2,".GseaPreranked(.*)"),full.names=T)
  result_file <- list.files(path=result_dir,pattern="gsea_report_for_na_pos_(.*).xls",full.names=T)
  gsea_result <- read.table(result_file,header = T,sep="\t",row.names=1)
  gsea_pathways <- str_to_title(rownames(gsea_result))
  gsea_pathways <- str_replace(gsea_pathways,"Tca","TCA")
  gsea_pathways <- str_replace(gsea_pathways,"Gpi","GPI")
  enrich_data_df <- rbind(enrich_data_df,data.frame(x=t2,y=gsea_pathways,NES=gsea_result$NES,PVAL=gsea_result$NOM.p.val))
}

#remove pvalue <0.01 pathways
min_pval <- by(enrich_data_df$PVAL, enrich_data_df$y, FUN=min)
select_pathways <- names(min_pval)[(min_pval<=0.01)]
select_enrich_data_df <- enrich_data_df[enrich_data_df$y%in% select_pathways,]
#converto pvalue to -log10
pvals <- select_enrich_data_df$PVAL
pvals[pvals<=0] = 1e-10
select_enrich_data_df$PVAL <- -log10(pvals)

#sort
pathway_pv_sum <- by(select_enrich_data_df$PVAL,select_enrich_data_df$y,FUN=sum)
pathway_order <- names(pathway_pv_sum)[order(pathway_pv_sum,decreasing = T)]
###########################top 10
##check before doing this 
pathway_order <- pathway_order[1:10]
select_enrich_data_df <- select_enrich_data_df[select_enrich_data_df$y %in% pathway_order,]
########################################
select_enrich_data_df$y <- factor(select_enrich_data_df$y,levels = pathway_order)

# #buble plot
p <- ggplot(select_enrich_data_df, aes(x = x, y = y, size = PVAL, color = NES)) +
  geom_point(shape=19) +
  #ggtitle("pathway heterogeneity") +
  labs(x = NULL, y = NULL,
       size = "-log10 pvalue", color = "NES") +
  scale_size(range = c(0, 2.5)) +
  scale_color_gradient( low = "white", high = "red") +
  #scale_color_gradient2(low="red",mid="white",high="blue",midpoint = 1) +
  theme(legend.position = "bottom", legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.key.size = unit(0.1, "cm"),
        legend.text = element_text(colour="black",size=6),
        axis.line = element_line(size=0.3, colour = "black"),
        #panel.grid.major = element_line(colour = "#d3d3d3"),
        #panel.grid.minor = element_blank(),
        axis.ticks = element_line(colour = "black", size = 0.3),
        panel.border = element_blank(), panel.background = element_blank(),
        axis.text.x=element_text(colour="black", size = 6,angle=90,hjust=1,vjust=0.5),
        axis.text.y=element_text(colour="black", size = 6)) +
  theme(plot.margin = unit(rep(1,4),"lines"))
ggsave(file.path(outDir,"non-malignant_enriched_pathway.pdf"),p,width = 3.6,height=2.5,units="in",device="pdf",useDingbats=FALSE)

##plot variance
p <- ggplot(pc_plotdata) + geom_point(aes(x,y,colour=factor(sel)),size=0.5) +
  scale_color_manual(values=c("gray","#ff4000")) +
  facet_wrap(~factor(types),scales="free",ncol = 4) + theme_bw() + 
  labs(x="Principal components", y="Explained variance (%)") +
  theme(legend.position="none",panel.grid.major = element_blank(), 
        panel.grid.minor= element_blank(),
        axis.line=element_line(size=0.2,colour="black"),
        axis.ticks = element_line(colour = "black",size=0.2),
        axis.text.x=element_text(colour="black", size = 6),
        axis.text.y=element_text(colour="black", size = 6),
        strip.background = element_rect(fill="white",size=0.2,colour = NULL),
        strip.text=element_text(size=6))

ggsave(file.path(outDir,"non-malignant_PC_variance_plot.pdf"),p,width = 7.5,height=2.7,units="in",device="pdf",useDingbats=FALSE)
unlink("preRankResults",recursive=T)
unlink("prerank.rnk")
date_string <- Sys.Date()
date_split <- strsplit(as.character(date_string),"-")[[1]]
unlink(paste0(tolower(month.abb[as.numeric(date_split[2])]),date_split[3]),recursive=T)