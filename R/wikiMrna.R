#' @title wikiMrna
#' @description Identify genes that are in common in both the wikipathway of
#' interest and the significantly differentially expressed input mRNAs.
#' @param MAE MultiAssayExperiment which will store the results of wikiMrna.
#' It is recommended to use the same MAE which stores output from the
#' diffExpressRes and reduceWiki functions.
#' @param mRNA_express Dataframe from the diffExpressRes function used on the
#' input mRNA data. This should be found as an assay within the MAE used in the
#' diffExpressRes function.
#' @param singleWiki Dataframe containing information about only one pathway.
#' This is output from the reduceWiki function. This should be found as an
#' assay within the MAE used in the reduceWiki function.
#' @param stringWiki Name of the pathway of interest. Should be the same
#' as the stringWiki parameter from the reduceWiki function.
#' @return A dataframe which only contains mRNAs which are found in both the
#' input data and the wikipathway of interest. Output will be stored as an assay
#' in the input MAE.
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
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus', org.Mm.eg.db)
#'
#' MAE <- dloadGmt(MAE = MAE, species = "Mus musculus")
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

    if (missing(MAE)) stop("MAE is missing. Add MAE. This will store the results from wikiMrna. Please use diffExpressRes and reduceWiki first.")

    if (missing(mRNA_express)) stop("mRNA_express is missing. Add dataframe which contains mRNA results from DE (e.g. log2fc, aveExpression) and gene IDs. Please use the diffExpressRes function on input mRNA data first. Output of diffExpressRes should be stored as an assay within the MAE used in the diffExpressRes function.")

    if (missing(singleWiki)) stop("singleWiki is missing. Add Dataframe which has data on which mRNAs are found within a single wikipathway. Please use the reduceWiki function first. Output of reduceWiki should be stored as an assay in the MAE used in the reduceWiki function.")

    if (missing(stringWiki)) stop("stringWiki is missing. Add Name of the wikipathway. Should be same as the stringWiki used in reduceWiki.")

    # Take dataframes out of MAE objects
    mRNAs <- mRNA_express

    pathway <- singleWiki

    # Find which genes from data and from selected pathway are the same
    GenesofInterest <- mRNAs[which(mRNAs$ID %in% pathway[,2]),]

    MAE2 <- suppressWarnings(suppressMessages(MultiAssayExperiment(list(x = GenesofInterest))))

    # unique name for each pathway
    a <- "GoI_"

    names(MAE2) <- paste0(a, stringWiki)

    MAE <- suppressWarnings(c(MAE, MAE2))

    return(MAE)
}
