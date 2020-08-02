#' @title combineGenes
#' @description Combines miR and mRNA data into one dataframe. Input columns
#' should be written as :timepoint.DifferentialExpressionResultType e.g.
#' D1.log2fc or H6.adjPval. Column names should be the same for miR and mRNA
#' data. A more detailed explanation of column nomenclature can be found in the
#' vignette. Resulting information will be stored as an assay within a MAE
#' object. combineGenes is essential for combined analysis of miR-mRNA data.
#' If using separate analysis, there is no need to use combineGenes.
#' @param MAE Input MAE object will have the results of combineGenes added to
#' it as an assay. It is recommended to use the MAE object which stores output
#' from startObject.
#' @param miR_data microRNA dataframe stored as an assay in a MAE. Rows
#' should be genes, columns are DE results and time point. This should be
#' the stored within the MAE used in the startObject function.
#' @param mRNA_data mRNA dataframe stored as an assay in a MAE. Rows
#' should be genes, columns are DE results and time point. This should be
#' the stored within the MAE used in the startObject function.
#' @return A dataframe which combines miR and mRNA data.
#' @export
#' @importFrom gtools mixedsort
#' @usage combineGenes(MAE, miR_data, mRNA_data)
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- combineGenes(MAE = MAE, miR_data = assay(MAE, 1),
#'                     mRNA_data = assay(MAE, 2))
combineGenes <- function(MAE, miR_data, mRNA_data){

    if (missing(MAE)) stop('Add MultiAssayExperiment to have the output of
                           combineGenes stored. Please use the
                           startObject function first.')

    if (missing(miR_data)) stop('Add dataframe of miR data. Rows should be
                                genes, columns are DE results and time point.
                                Please use the startObject function first.
                                Output of the startObject function will be
                                stored as assays within the MAE used in the
                                startObject function.')

    if (missing(mRNA_data)) stop('Add dataframe of mRNA data. Rows should be
                                genes, columns are DE results and time point.
                                Please use the startObject function first.
                                Output of the startObject function will be
                                stored as assays within the MAE used in the
                                startObject function.')

    # Extract the miR and mRNA data
    miR_data <- as.data.frame(miR_data)

    mRNA_data <- as.data.frame(mRNA_data)

    # Order them based on time point and result type
    miR_order <- gtools::mixedsort(names(miR_data))

    mRNA_order <- gtools::mixedsort(names(mRNA_data))

    # New order is how both data sets should be ordered
    miR_data <- miR_data[miR_order]

    mRNA_data <- mRNA_data[mRNA_order]

    genetic_data <- as.data.frame(rbind(miR_data, mRNA_data))

    # Combine with MAE object
    MAE <- c(MAE, suppressMessages(MultiAssayExperiment(experiments = list(
                                              "genetic_data" = genetic_data))))

    return(MAE)
}
