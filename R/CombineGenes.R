#' @title CombineGenes
#' @description Combines miR and mRNA data into one dataframe. Input columns
#' should be written as follows:
#' timepoint.DEresult e.g. D1.log2fc or H6.adjPval.
#'
#' @param miR_data microRNA dataframe
#' @param mRNA_data mRNA dataframe
#'
#' @return A dataframe with combines miR and mRNA data.
#' @export
#' @import gtools
#' @usage CombineGenes(miR_data, mRNA_data)
#' @examples
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' CombineGenes(miR_data = miR, mRNA_data = mRNA) -> genetic_data
CombineGenes <- function(miR_data, mRNA_data){
        gtools::mixedsort(names(miR_data)) -> miR_order
        gtools::mixedsort(names(mRNA_data)) -> mRNA_order
        miR_data <- miR_data[miR_order]
        mRNA_data <- mRNA_data[mRNA_order]
        rbind(miR_data, mRNA_data) -> genetic_data
        return(genetic_data)
}
