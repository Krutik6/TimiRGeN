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
#' StartObject(miR = NULL, mRNA = mRNA) -> MAE
#' getIDs_mRNA_mouse(MAE, MAE@ExperimentList$mRNA) -> MAE
#' Express(df = MAE@ExperimentList$mRNA, dataType = 'Log2FC', 
#' genes_ID = MAE@ExperimentList$mRNA_entrez,
#' idColumn = 'GENENAME') -> mRNA_express
# 
#' dloadGMT(MAE = MAE, speciesInitials = "Mm") -> MAE
# 
#' ReduceWiki(path_data = MAE@ExperimentList$path_data,
#' stringWiki = 'TGF Beta Signaling Pathway') -> singlewiki
# 
#' Express(df = MAE@ExperimentList$mRNA, dataType = 'Log2FC',
#' genes_ID = MAE@ExperimentList$mRNA_entrez,
#' idColumn = 'GENENAME') -> MAE@ExperimentList$mRNA_express
# 
#' WikimRNA(mRNA_express = MAE@ExperimentList$mRNA_express,
#' SingleWiki = singlewiki) -> MAE@ExperimentList$GenesofInterest
WikimRNA <- function(mRNA_express, SingleWiki){
if (missing(mRNA_express)) stop("Use Express function of mRNA data.")
if (missing(SingleWiki)) stop("Use ReduceWiki function on a Wikipathway
of interest.")
mRNA_express[which(mRNA_express$ID %in% SingleWiki[,2]),] -> GenesofInterest
return(GenesofInterest)
}
