#' @title startObject
#' @description Creates a MultiAssayExperiment (MAE) from miR and mRNA
#' dataframes. MAE's will be the constant object used throughout TimiRGeN. The
#' input dataframes should contain rows as genes, and results from differential
#' expression (DE) as columns. Columns should also indicate the time point
#' related to each sample. Row names and column names must adhere to TimiRGeN
#' friendly nomenclature. Please do read the vignette for a full description of
#' the required nomenclature.
#' @param miR microRNA dataframe/ matrix. Rows should be miR gene names which
#' use the TimiRGeN friendly naming system. Columns should be results of
#' DE and time points.
#' @param mRNA mRNA dataframe/ matrix. Rows should be mRNA gene names.
#' Columns should be results of DE and time points.
#' @return MultiAssayExperiment containing miR and mRNA data stored as assays.
#' Use assays(MAE, i) or MAE[[i]] to access assays. Use metadata(MAE)[[i]] to
#' access metadata. Use experiments(MAE)[[i]] to access experiments.
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

    if (missing(miR)) stop('
                           miR is missing.
                           Add microRNAs. Rows should be gene names and columns
                           should be DE results. Nomenclature is explained in
                           detail within the vignette.')

    if (missing(mRNA)) stop('
                            mRNA is missing.
                            Add mRNAs. Rows should be gene names and columns
                            should be DE results. Nomenclature is explained in
                            detail within the vignette.')

    Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))

    MAE <- suppressMessages(MultiAssayExperiment(experiments = Data))

    return(MAE)
}
