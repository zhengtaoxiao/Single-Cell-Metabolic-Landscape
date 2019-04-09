##gmtPathways is copied from fgsea package.
gmtPathways <- function(gmt.file) {
    pathwayLines <- strsplit(readLines(gmt.file), "\t")
    pathways <- lapply(pathwayLines, tail, -2)
    names(pathways) <- sapply(pathwayLines, head, 1)
    pathways
}

## get_os is copied from https://github.com/r-lib/rappdirs/blob/master/R/utils.r#L1
get_os <- function() {
  if (.Platform$OS.type == "windows") { 
    "win"
  } else if (Sys.info()["sysname"] == "Darwin") {
    "mac" 
  } else if (.Platform$OS.type == "unix") { 
    "unix"
  } else {
    stop("Unknown OS")
  }
}

##from :  https://github.com/mikelove/DESeq2/blob/master/R/core.R
estimateSizeFactorsForMatrix <- function(counts, locfunc=stats::median,
                                         geoMeans, controlGenes) {
  if (missing(geoMeans)) {
    incomingGeoMeans <- FALSE
    loggeomeans <- rowMeans(log(counts))
  } else {
    incomingGeoMeans <- TRUE
    if (length(geoMeans) != nrow(counts)) {
      stop("geoMeans should be as long as the number of rows of counts")
    }
    loggeomeans <- log(geoMeans)
  }
  if (all(is.infinite(loggeomeans))) {
    stop("every gene contains at least one zero, cannot compute log geometric means")
  }
  sf <- if (missing(controlGenes)) {
    apply(counts, 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeans)[is.finite(loggeomeans) & cnts > 0]))
    })
  } else {
    if ( !( is.numeric(controlGenes) | is.logical(controlGenes) ) ) {
      stop("controlGenes should be either a numeric or logical vector")
    }
    loggeomeansSub <- loggeomeans[controlGenes]
    apply(counts[controlGenes,,drop=FALSE], 2, function(cnts) {
      exp(locfunc((log(cnts) - loggeomeansSub)[is.finite(loggeomeansSub) & cnts > 0]))
    })
  }
  if (incomingGeoMeans) {
    # stabilize size factors to have geometric mean of 1
    sf <- sf/exp(mean(log(sf)))
  }
  sf
}


##calculate how many pathways of one gene involved.
num_of_pathways <- function (gmtfile,overlapgenes){
  pathways <- gmtPathways(gmtfile)
  pathway_names <- names(pathways)
  filter_pathways <- list()
  for (p in pathway_names){
    genes <- pathways[[p]]
    common_genes <- intersect(genes,overlapgenes)
    if(length(common_genes>=5)){
      filter_pathways[[p]] <- common_genes
    }
  }
  
  all_genes <- unique(as.vector(unlist(filter_pathways)))
  gene_times <- data.frame(num =rep(0,length(all_genes)),row.names = all_genes)
  for(p in pathway_names){
    for(g in filter_pathways[[p]]){
      gene_times[g,"num"] = gene_times[g,"num"]+1
    }
  }
  gene_times
} 
