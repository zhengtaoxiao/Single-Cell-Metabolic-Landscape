library(scater)
library(stringr)
options(stringsAsFactors=FALSE)
library(RColorBrewer)
source("../utils.R")

pathway_file <- "../Data/KEGG_metabolism.gmt"
hallmark_gmt <- '../Data/h.all.v6.1.symbols.gmt'

pathways <- gmtPathways(pathway_file)
all_pathways <- gmtPathways(hallmark_gmt)

oxphos_genes <- pathways[["Oxidative phosphorylation"]]
glycolysis_genes <- pathways[["Glycolysis / Gluconeogenesis"]]
hypoxia_genes <- all_pathways[["HALLMARK_HYPOXIA"]]

args <- commandArgs()
tumor <- args[6]
outDir <- file.path("dataset",tumor)
if(!dir.exists(outDir) ) dir.create(outDir,recursive=TRUE)

#1. Loading the data
selected_sce <- readRDS(file.path("../1-ReadData/dataset/",tumor,"selected_sce.rds"))
selected_tumor_sce <- selected_sce[,selected_sce$cellType=="Malignant"]
#=========================================================================
tumors <- unique(selected_tumor_sce$tumor)
gene_pathway_num <- num_of_pathways(pathway_file,intersect(unlist(pathways),rownames(selected_sce)))

all_exp <- assay(selected_tumor_sce,"exprs")

oxphos_exp <- all_exp[rownames(all_exp)%in% oxphos_genes,]
glycolysis_exp <- all_exp[rownames(all_exp)%in% glycolysis_genes,]
hypoxia_exp <- all_exp[rownames(all_exp)%in% hypoxia_genes,]

oxphos <- colMeans(as.matrix(oxphos_exp),na.rm=T)
glycolysis <- colMeans(as.matrix(glycolysis_exp),na.rm=T)
hypoxia <- colMeans(as.matrix(hypoxia_exp),na.rm=T)
dat <- data.frame(OXPHOS=oxphos,Glycolysis=glycolysis,Hypoxia=hypoxia)

print("correlation:")
print(cor(dat))

#correlation plot for each two of them
dat_min <- 0
dat_max <- 4
p=ggplot(dat,aes(x=OXPHOS,y=Glycolysis)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  labs(x = "OXPHOS", y = "Glycolysis") +
  theme(axis.line=element_line(size=0.3,colour="black"),
       axis.ticks = element_line(size=0.3,color="black"),
       axis.text.x=element_text(size=6),
       axis.text.y=element_text(size=6),
       axis.title.x=element_text(size=8),
       axis.title.y=element_text(size=8))

ggsave(filename = file.path(outDir,"malignant_oxphos_glycolysis.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

p=ggplot(dat,aes(x=OXPHOS,y=Hypoxia)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  labs(x = "OXPHOS", y = "Hypoxia") +
  theme(axis.line=element_line(size=0.3,colour="black"),
        axis.ticks = element_line(size=0.3,color="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))
  

ggsave(filename = file.path(outDir,"malignant_oxphos_hypoxia.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

p=ggplot(dat,aes(x=Glycolysis,y=Hypoxia)) + 
  geom_point(size=0.5) +
  geom_smooth(method="lm",color="red") +
  xlim(dat_min,dat_max) + ylim(dat_min,dat_max) +
  labs(x = "Glycolysis", y = "Hypoxia") +
  theme_classic()  + theme(aspect.ratio = 0.8) +
  theme(axis.line=element_line(size=0.3,colour="black"),
        axis.ticks = element_line(size=0.3,color="black"),
        axis.text.x=element_text(size=6),
        axis.text.y=element_text(size=6),
        axis.title.x=element_text(size=8),
        axis.title.y=element_text(size=8))

ggsave(filename = file.path(outDir,"malignant_glycolysis_hypoxia.pdf"),p,device = "pdf",width=2,height=1.5,units="in",useDingbats=FALSE)

# #correlation in each tumor
# for( i in tumors){
#   each_exp <- assay(selected_tumor_sce[,selected_tumor_sce$tumor==i],"exprs")
#   oxphos_exp <- each_exp[rownames(each_exp)%in% oxphos_genes,]
#   glycolysis_exp <- each_exp[rownames(each_exp)%in% glycolysis_genes,]
#   hypoxia_exp <- each_exp[rownames(each_exp)%in% hypoxia_genes,]
  
#   oxphos <- colMeans(as.matrix(oxphos_exp),na.rm=T)
#   glycolysis <- colMeans(as.matrix(glycolysis_exp),na.rm=T)
#   hypoxia <- colMeans(as.matrix(hypoxia_exp),na.rm=T)
#   dat <- data.frame(OXPHOS=oxphos,Glycolysis=glycolysis,Hypoxia=hypoxia)
#   print(paste0("correlations in patient:",i))
#   print(cor(dat))
# }


#non-tumors 
selected_nontumor_sce <- selected_sce[,selected_sce$cellType!="Malignant"]
cell_types <- unique(selected_nontumor_sce$cellType)
cor_matrix <- matrix(NA,nrow=length(cell_types),ncol=3,dimnames = list(cell_types,c("oxphos_glycolosis","oxphos_hypoxia","glycolysis_hypoxia")))
for(c in cell_types){
  each_exp <- assay(selected_nontumor_sce[,selected_nontumor_sce$cellType==c],"exprs")
  oxphos_exp <- each_exp[rownames(each_exp)%in% oxphos_genes,]
  glycolysis_exp <- each_exp[rownames(each_exp)%in% glycolysis_genes,]
  hypoxia_exp <- each_exp[rownames(each_exp)%in% hypoxia_genes,]
  
  oxphos <- colMeans(as.matrix(oxphos_exp),na.rm=T)
  glycolysis <- colMeans(as.matrix(glycolysis_exp),na.rm=T)
  hypoxia <- colMeans(as.matrix(hypoxia_exp),na.rm=T)
  dat <- data.frame(OXPHOS=oxphos,Glycolysis=glycolysis,Hypoxia=hypoxia)
  #calculate correlation
  cor_matrix[c,1] <- cor(dat[,1:2])[1,2]
  cor_matrix[c,2] <- cor(dat[,c(1,3)])[1,2]
  cor_matrix[c,3] <- cor(dat[,2:3])[1,2]
}
write.table(cor_matrix,file.path(outDir,"non-malignant_Oxphos_Glycolysis_Hypoxia_cor.txt"),row.names = T,col.names=T,sep="\t")
