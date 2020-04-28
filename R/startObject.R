#' @title startObject
#' @aliases startObject
#' @description Creates a MultiAssayObject from miR and mRNA dataframes. This
#' will be the constant object used throughout TimiRGeN.
#' @param miR microRNA dataframe/ matrix
#' @param mRNA mRNA dataframe/ matrix
#' @return MultiAssayExperiment object
#' @export
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @usage startObject(miR, mRNA)
#' @examples
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' MAE <- startObject(miR = miR, mRNA = mRNA)
startObject <- function(miR, mRNA){
    Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))
    MAE <- suppressMessages(MultiAssayExperiment(experiments = Data))
    return(MAE)
}
