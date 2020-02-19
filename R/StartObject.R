#' @title StartObject
#' @aliases StartObject
#' @description Creates a MultiAssayObject from miR and mRNA dataframes. This
#' will be the constant object used throughout TimiRGeN.
#' @param miR microRNA dataframe/ matrix
#' @param mRNA mRNA dataframe/ matrix
#'
#' @return MultiAssayExperiment object
#' @export 
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @usage StartObject(miR, mRNA)
#' @examples
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
StartObject <- function(miR, mRNA){
Data <- list("miR" = as.data.frame(miR), "mRNA" = as.data.frame(mRNA))
MAE <- MultiAssayExperiment(experiments = Data)
return(MAE)
}