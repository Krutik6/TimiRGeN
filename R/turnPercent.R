#' @title turnPercent
#' @description Genes found in common between the input data and each
#' pathway are normalised by percentages. This is to normalise for pathway size.
#' @param MAE MultiAssayExperiment which will store the output from
#' turnPercent. It is recommended to use the MAE which stores output from the
#' wikiMatrix function.
#' @param wikiMatrix Numeric matrix of wikipathways and samples. This should be
#' stored as an assay within the MAE used in the wikiMatrix function.
#' @return A percentage matrix which contrasts genes found in pathways and
#' samples. Output will be stored as an assay within the input MAE.
#' @export
#' @usage turnPercent(MAE, wikiMatrix)
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
#' MAE <- turnPercent(MAE = MAE, wikiMatrix = assay(MAE, 1))
turnPercent <- function(MAE, wikiMatrix){

    if (missing(MAE)) stop('MAE is missing. Add MAE. Output of turnPercent will be stored within this MAE. Please use wikiMatrix first.')

    if (missing(wikiMatrix)) stop('wikiMatrix is missing. Add matrix of wikipathways and samples. Please use the wikiMatrix function before using turnPercent. Output of wikiMatrix should be stored as an assay within the MAE used in the wikiMatrix function.')

    rowInt <- nrow(wikiMatrix)

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
