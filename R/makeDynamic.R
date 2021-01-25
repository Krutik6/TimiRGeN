#' @title makeDynamic
#' @description Produces a dataframe that can be imported into pathvisio
#' to show how changes in genes expression levels over the time course. Follow
#' instructions found in the vignette which show on how to save this file and
#' further instructions found in /inst/Pathvisio_GRN_guide.pdf to see how
#' this can help in GRN construction.
#' @param MAE MultiAssayExpreriment to store the output of makeDynamic. It is
#' recommended to use the same MAE which stores output from matrixFilter.
#' @param miR_expression Dataframe containing abundance values
#' (e.g. log2fc or average expression) from miR specific differential expression
#' , along with gene IDs. This is the output from diffExpressRes. Output of
#' diffExpressRes should be stored as an assay within the MAE used in
#' diffExpressRes.
#' @param mRNA_expression Dataframe containing abundance values (log2fc or
#' average expression) from mRNA specific differential expression, along with
#' gene IDs. This is the output from diffExpressRes. Output of diffExpressRes
#' should be stored as an assay within the MAE used in diffExpressRes.
#' @param miR_IDs_adj Dataframe which contain adjusted gene IDs from miR data.
#' Either miR_adjusted_entrez or miR_adjusted_ensembl. Should be found
#' as an assay in the MAE used a getIdsMir function.
#' @param dataType String which represents the gene ID used in this analysis.
#' Either "En" (ensembl data) or "L" (entrez data).
#' @return miR and mRNA dynamic data that can be saved and be used in
#' pathvisio to display dynamic behaviour of miRs and mRNAs of interest over
#' the time series in a signalling pathway of interest. Output will be stored
#' as an assay in the input MAE.
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
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus')
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

    if (missing(MAE)) stop('
                           MAE is missing.
                           Add MAE. This will store the output of
                           makeDynamic. Please use matrixFilter first.')

    if (missing(miR_expression)) stop('
                                      miR_expression is missing.
                                      Add dataframe which contains DE
                                      abundance values and gene IDs. Please
                                      use diffExpressRes on miR data first.
                                      Output of diffExpressRes should be stored
                                      as an assay within the MAE used in the
                                      diffExpressRes function.')

    if (missing(mRNA_expression)) stop('
                                      mRNA_expression is missing.
                                      Add dataframe which contains DE
                                      abundance values and gene IDs. Please
                                      use diffExpressRes on mRNA data first.
                                      Output of diffExpressRes should be stored
                                      as an assay within the MAE used in the
                                      diffExpressRes function.')

    if (missing(miR_IDs_adj)) stop('
                                   miR_IDs_adj is missing.
                                   Add adjusted miR gene ID data.
                                   Please use the getIdsMirHuman or
                                   getIdsMirMouse function first. Output of
                                   a getIdsMir function should be stored as
                                   assays within the MAE used in the
                                   getIds function.')

    if (missing(dataType)) stop('
                                dataType is missing.
                                Add a string. "En" for ensembl or "L"
                                for entrez.')

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
