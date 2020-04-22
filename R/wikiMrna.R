#' @title wikiMrna
#' @description Identify which mRNAs are in common with the genes found in the
#' wikipathway of interest.
#' @param MAE MultiAssayExperiment object.
#' @param mRNA_express Output from the Express function on mRNA data.
#' @param singleWiki Output from ReduceWiki function on a Wikipathway
#' of interest.
#' @param stringWiki Name of the pathway of interest.
#' @return A dataframe of mRNA entries which are also found in the wikipathway
#' of interest.
#' @export
#' @usage wikiMrna(MAE, mRNA_express, singleWiki, stringWiki='')
#' @examples
#' library(biomaRt)
#' library(MultiAssayExperiment)
#' miR <- mm_miR[1:100,]
#' mRNA <- mm_mRNA[1:200,]
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
#' MAE <- getIdsMrnaMouse(MAE, assay(MAE, 2), "useast")
#' MAE <- dloadGmt(MAE = MAE, speciesInitials = "Mm")
#'
#' MAE <- reduceWiki(MAE, path_data = assay(MAE, 11),
#'                   stringWiki = 'TGF Beta Signaling Pathway')
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 7),
#'                idColumn = 'GENENAME',
#'                name = "mRNA_log")
#'
#' MAE <- wikiMrna(MAE, mRNA_express = assay(MAE, 13),
#'                 singleWiki = assay(MAE, 12),
#'                 stringWiki = 'TGF Beta Signaling Pathway')
wikiMrna <- function(MAE, mRNA_express, singleWiki, stringWiki){

    if (missing(MAE)) stop("Insert MAE object.");
    if (missing(mRNA_express)) stop("Use Express function of mRNA data.");
    if (missing(singleWiki)) stop("Use ReduceWiki function on a Wikipathway
                                   of interest.");
    if (missing(stringWiki)) stop("Name of the wikipathway.");

    mRNAs <- mRNA_express
    pathway <- singleWiki
    GenesofInterest <- mRNAs[which(mRNAs$ID %in% pathway[,2]),]
    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = GenesofInterest)))
    a <- "GoI_"
    names(MAE2) <- paste0(a, stringWiki)
    MAE <- c(MAE, MAE2)
    return(MAE)
}
