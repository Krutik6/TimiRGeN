#' @title combineGenes
#' @description Combines miR and mRNA data into one dataframe. Input columns
#' should be written as :timepoint.Differentialexpressionresulttyle e.g.
#' D1.log2fc or H6.adjPval.
#' @param MAE MultiAssayExperiment
#' @param miR_data microRNA dataframe
#' @param mRNA_data mRNA dataframe
#' @return A dataframe with combines miR and mRNA data which can be stored
#' within the MAE object.
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

    if (missing(MAE)) stop('Add MultiAssayExperiment.')

    if (missing(miR_data)) stop('miR data from startObject.')

    if (missing(mRNA_data)) stop('mRNA data from startObject.')

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
