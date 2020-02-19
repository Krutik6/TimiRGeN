#' @title CombineGenes
#' @description Combines miR and mRNA data into one dataframe. Input columns
#' should be written as follows:
#' timepoint.DEresult e.g. D1.log2fc or H6.adjPval.
#'
#' @param miR_data microRNA dataframe
#' @param mRNA_data mRNA dataframe
#'
#' @return A dataframe with combines miR and mRNA data which can be stored 
#' within the MAE object.
#' @export
#' @import gtools
#' @usage CombineGenes(miR_data, mRNA_data)
#' @examples
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' CombineGenes(miR_data = MAE@ExperimentList$miR, mRNA_data = 
#' MAE@ExperimentList$mRNA) -> MAE@ExperimentList$genetic_data
CombineGenes <- function(miR_data, mRNA_data){
miR_data <- as.data.frame(miR_data)
mRNA_data <- as.data.frame(mRNA_data)
miR_order <- gtools::mixedsort(names(miR_data))
mRNA_order <- gtools::mixedsort(names(mRNA_data)) 
miR_data <- miR_data[miR_order]
mRNA_data <- mRNA_data[mRNA_order]
genetic_data <- as.data.frame(rbind(miR_data, mRNA_data))
return(genetic_data)
}
