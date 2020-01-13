#' @title WikimRNA
#' @description Identify which mRNAs are in common with the genes found in the
#' wikipathway of interest.
#' @param mRNA_express Output from the Express function on mRNA data.
#' @param SingleWiki Output from ReduceWiki function on a Wikipathway
#' of interest.
#' @return A dataframe of mRNA entries which are also found in the wikipathway
#' of interest.
#' @export
#' @usage WikimRNA(mRNA_express, SingleWiki)
#' @examples
#' library(biomaRt)
#' mm_mRNA -> mRNA
#' getIDs_mRNA_mouse(mRNA)
#' Express(df = mRNA, dataType = 'Log2FC', genes_ID = mRNA_entrez,
#' idColumn = 'GENENAME') -> mRNA_express
#' path_data <- data.frame("wpid" = c(rep("WP571", 6)),
#' "gene" = c(16175, 12370,26419, 19249, 19645, 18479),
#' "name" = c(rep("Fas pathway and Stress induction of HSP regulation", 6)))
#' ReduceWiki(path_data = path_data,
#' stringWiki = 'Fas pathway and Stress induction of
#' HSP regulation') -> singlewiki
#' WikimRNA(mRNA_express = mRNA_express,
#' SingleWiki = singlewiki) -> GenesofInteres
WikimRNA <- function(mRNA_express, SingleWiki){
if (missing(mRNA_express)) stop("Use Express function of mRNA data.")
if (missing(SingleWiki)) stop("Use ReduceWiki function on a Wikipathway
of interest.")
mRNA_express[which(mRNA_express$ID %in% SingleWiki[,2]),] -> GenesofInterest
return(GenesofInterest)
}
