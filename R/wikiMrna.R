#' @title wikiMrna
#' @description Identify which genes are in common with the genes found in the
#' wikipathway of interest and the significantly differentially
#' expressed input mRNAs. Store the results in a MAE as an assay.
#' @param MAE MultiAssayExperiment object which will have the results of
#' wikiMrna added to it. It is recommended to use the same MAE which stores
#' information from diffExpressRes and reduceWiki functions.
#' @param mRNA_express Dataframe from the diffExpressRes function used on the
#' input mRNA data. This should be found as an assay within the MAE used in
#' diffExpressRes.
#' @param singleWiki Output from reduceWiki function on a Wikipathway
#' of interest. This information should be found as an assay within the MAE
#' used in the reduceWiki function.
#' @param stringWiki Name of the pathway of interest. Should be the same
#' as the stringWiki from the reduceWiki function.
#' @return A dataframe which only contains mRNAs from the input data and
#' the wikipathway of interest.
#' @export
#' @usage wikiMrna(MAE, mRNA_express, singleWiki, stringWiki='')
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR[1:200,]
#'
#' mRNA <- mm_mRNA[401:600,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
#'
#' MAE <- getIdsMrnaMouse(MAE, assay(MAE, 2), "www")
#'
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

    if (missing(MAE)) stop("Add MAE. This will store the results
                           from wikiMrna in it. Please use the
                           diffExpressRes and reduceWiki first.")

    if (missing(mRNA_express)) stop("Add dataframe which contains mRNAs,
                                    results from DE e.g. log2fc, aveExpression,
                                    and gene IDs. Please use the diffExpressRes
                                    function on input mRNA data first.
                                    Output of diffExpressRes should be stored
                                    as an assay in the MAE used in the
                                    diffExpressRes function.")

    if (missing(singleWiki)) stop("Add Dataframe which has data about which
                                  mRNAs are found within a single wikipathway.
                                  Please use the reduceWiki function first.
                                  Output of singeWiki should be stored as an
                                  assay in the MAE used in the singleWiki
                                  function.")

    if (missing(stringWiki)) stop("Add Name of the wikipathway. Should be same
                                  as the stringWiki used in reduceWiki.")

    # Take dataframes out of MAE objects
    mRNAs <- mRNA_express

    pathway <- singleWiki

    # Find which genes from data and from selected pathway are the same
    GenesofInterest <- mRNAs[which(mRNAs$ID %in% pathway[,2]),]

    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = GenesofInterest)))

    # unique name for each pathway
    a <- "GoI_"

    names(MAE2) <- paste0(a, stringWiki)

    MAE <- c(MAE, MAE2)

    return(MAE)
}
