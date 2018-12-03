#install scater,scran,Rstne
source("https://bioconductor.org/biocLite.R")
biocLite("scater")
biocLite("scran")
biocLite("Rtsne")
biocLite("biomaRt")

#scimpute
if(library(devtools,logical.return=T)){
  install_github("Vivianstats/scImpute",force=TRUE)
}else{
  warning("please install devtools and try again")
  #install.packages("devtools")
}

#pheatmap
if(!library(pheatmap,logical.return=T)){
  warning("please install pheatmap and try again")
  #install.packages("pheatmap")
}

#ggrepel
if(!library(ggrepel,logical.return=T)){
  warning("please install pheatmap and try again")
  #install.packages("ggrepel")
}

# download gsea
if(!file.exists("gsea-3.0.jar")){
  warning("please download gsea-3.0.jar from http://software.broadinstitute.org/gsea/downloads.jsp")
}