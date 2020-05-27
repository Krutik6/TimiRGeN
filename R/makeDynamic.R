#' @title makeDynamic
#' @description Produces a data frame that can be imported into pathvisio
#' to show changes in genes over time. Follow vignette instructions on how to
#' save this file and further instructions found in the
#' /inst/Pathvisio_GRN_guide.pdf to see how this can help in GRN construction.
#' @param MAE MultiAssayObject.
#' @param miR_expression microRNA data after going through diffExpressRes.
#' @param mRNA_expression mRNA data after going through diffExpressRes.
#' @param miR_IDs_adj Either miR_adjusted_entrez or miR_adjusted_ensembl.
#' @param dataType Either En (ensembl data) or L (entrez data).
#' @return MicroRNA and mRNA dynamic data that can be used in pathvisio if
#'exported.
#' @export
#' @usage makeDynamic(MAE, miR_expression, mRNA_expression, miR_IDs_adj,
#'                    dataType = '')
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR[1:100,]
#'
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
#'
#' MAE <- getIdsMrnaMouse(MAE, assay(MAE, 2), "useast")
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 1), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 3),
#'                idColumn = 'GENENAME',
#'                name = "miR_express")
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 7),
#'                idColumn = 'GENENAME',
#'                name = 'mRNA_express')
#'
#' MAE <- makeDynamic(MAE, miR_expression = assay(MAE, 9),
#'                   mRNA_expression = assay(MAE, 10),
#'                   miR_IDs_adj = assay(MAE, 5),
#'                   dataType = "L")
makeDynamic <- function(MAE, miR_expression, mRNA_expression, miR_IDs_adj,
                        dataType){

    if (missing(MAE)) stop('Add MAE object')

    if (missing(miR_expression)) stop('Input miR expression data from
                                       diffExpressRes function.')

    if (missing(mRNA_expression)) stop('Input mRNA expression data from
                                        diffExpressResfunction.')

    if (missing(miR_IDs_adj)) stop('Input miR ID data that is adjusted for
                                    repeats. Ensembl or entrez.')

    if (missing(dataType)) stop('En for ensembl or L for entrez.')

    # retrieve data frames from MAE objects
    miR_expression <- miR_expression

    mRNA_expression <- mRNA_expression

    miR_IDs_adj <- miR_IDs_adj

    miR_expression$names <- rownames(miR_expression)

    # Merge time series data with IDs
    X <- merge(x= miR_expression, y= miR_IDs_adj, by.x= 'names',
               by.y= 'GENENAME', all=TRUE)

    rownames(X) <- X$names

    X$names <- X$ID.x <- NULL

    # tidy up data frame
    colnames(X) <- gsub(colnames(X), pattern = "ID.y", replacement = "ID")

    names(mRNA_expression) <- names(X)

    Dynamic <- rbind(X, mRNA_expression)

    # Add L or En as a new column to data frame
    Dynamic <- cbind(Dynamic, dataType)

    MAE2 <- suppressMessages(MultiAssayExperiment(list('Dynamic' = Dynamic)))

    MAE <- c(MAE, MAE2)

return(MAE)
}
