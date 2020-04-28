#' @title turnPercent
#' @description Converts integers into percentages, within a matrix.
#' @param MAE MultiAssayExperiment object.
#' @param wikiMatrix Matrix of wikipathways and samples.
#' @param rowInt Which row contains the total number of genes per wikipathway?
#' @return A percentage matrix.
#' @export
#' @usage turnPercent(MAE, wikiMatrix, rowInt)
#' @examples
#' MAE <- MultiAssayExperiment()
#' metadata(MAE)["e_list"] <- e_list
#' metadata(MAE)["w_list"] <- w_list[1:10]
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
#'
#' MAE <- turnPercent(MAE = MAE,
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
turnPercent <- function(MAE, wikiMatrix, rowInt){

    if (missing(MAE)) stop('Add MAE object');
    if (missing(wikiMatrix)) stop('Add matrix of wikipathways and samples.');
    if (missing(rowInt)) stop('Add the total number of samples.');

    df1 <- as.matrix(wikiMatrix)

    X <- t(t(df1)/df1[rowInt,]*100)
    X <- format(round(X, 2), nsmall = 2)

    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                             "Percentmatrix" = as.data.frame(X))))
    MAE <- c(MAE, MAE2)
return(MAE)
}
