#' @title startObject
#' @description Creates a MultiAssayObject from miR and mRNA dataframes. This
#' will be the constant object used throughout TimiRGeN. The dataframes should
#' contain rows as genes, and results from differential expression (DE) as
#' columns. The column names should be gene names must adhere to
#' TimiRGeN standards for processing and analysis, please read the vignette
#' for a full description of the required nomenclature.
#' @param miR microRNA dataframe/ matrix. Rows should be miR gene names which
#' should use TimiRGeN standard nomenclature. Columns should be results of
#' DE and time points.
#' @param mRNA mRNA dataframe/ matrix. Rows should be mRNA gene names.
#' Columns should be results of DE and time points.
#' @return MultiAssayExperiment object
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @usage startObject(miR, mRNA)
#' @examples
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
startObject <- function(miR, mRNA){

    if (missing(miR)) stop('Add microRNA. Rows should be gene names
                           and columns should be DE results. Nomenclature is
                           explained in detail within the vignette.')

    if (missing(mRNA)) stop('Add mRNA. Rows should be gene names
                           and columns should be DE results. Nomenclature is
                           explained in detail within the vignette.')

    Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))

    MAE <- suppressMessages(MultiAssayExperiment(experiments = Data))

    return(MAE)
}
