#' @title turnPercent
#' @description Genes found in common between the input data and each
#' pathway are normalised by percentages to normalise for pathway size. Stores
#' output as an assay in a MAE object.
#' @param MAE MultiAssayExperiment object which will have the output of
#' turnPercent added as an assay. It is recommended to use the MAE produced by
#' the wikiMatrix function.
#' @param wikiMatrix Matrix of wikipathways and samples. This should be stored
#' as an assay in the MAE used in the wikiMatrix function.
#' @param rowInt Which row contains the total number of genes per wikipathway?
#' This will be 1+ the number of samples in your input data. For example our
#' combined miR-mRNA analysis in the example has rowInt = 6 because our example
#' has 5 time points.
#' @return A percentage matrix.
#' @export
#' @usage turnPercent(MAE, wikiMatrix, rowInt)
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list_mouse
#'
#' metadata(MAE)[["w_list"]] <- w_list_mouse[1:10]
#'
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
#'
#' MAE <- turnPercent(MAE = MAE,
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
turnPercent <- function(MAE, wikiMatrix, rowInt){

    if (missing(MAE)) stop('Add MAE. Output of turnPercent will be
                           stored within this MAE. Please use wikiMatrix first.
                           ')

    if (missing(wikiMatrix)) stop('Add matrix of wikipathways and samples.
                                  Please use the wikiMatrix function before
                                  using turnPercent. Output of wikiMatrix
                                  should be stored as an assay within the MAE
                                  used in the wikiMatrix function.')

    if (missing(rowInt)) stop('Add an integer which represents the row that
                              contains the total number of genes.')

    df1 <- as.matrix(wikiMatrix)

    # Convert # genes found in pathway to percentage
    X <- t(t(df1)/df1[rowInt,]*100)

    # Round to 2 decimal places
    X <- format(round(X, 2), nsmall = 2)

    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                             "Percentmatrix" = as.data.frame(X))))

    MAE <- c(MAE, MAE2)

return(MAE)
}
