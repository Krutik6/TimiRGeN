#' @title reduceWiki
#' @description Return gene names for a single wikipathway of interest.
#' @param MAE MultiAssayExperiment object.
#' @param path_data File with wikipathway - gene data.
#' @param stringWiki Full name of the wikipathway of interest.
#' @return A dataframe which only contains information about the wikipathway of
#'interest.
#' @export
#' @usage reduceWiki(MAE, path_data, stringWiki = '')
#' @examples
#' library(MultiAssayExperiment)
#' MAE <- MultiAssayExperiment(list(path_data = data.frame(
#'                         "wpid" = c(rep("WP571", 6)),
#'                         "gene" = c(16175, 12370,26419, 19249, 19645, 18479),
#'                         "name" = c(rep("Fas pathway and Stress induction of HSP regulation",
#'                         6)))))
#' MAE <- reduceWiki(MAE, path_data = assay(MAE, 1),
#'                  stringWiki = 'Fas pathway and Stress induction of HSP regulation')
reduceWiki <- function(MAE, path_data, stringWiki){

    if (missing(MAE)) stop('Add a MultiAssayExperiment object.');
    if (missing(path_data)) stop('Input wikipathways - gene data.');
    if (missing(stringWiki)) stop('Input name of chosen wikipathway e.g.
        TGF-beta Signaling Pathway')

    singlewiki <- path_data[which(path_data$name == stringWiki),]
    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = singlewiki)))
    names(MAE2) <- stringWiki
    MAE <- suppressMessages(c(MAE, MAE2))
    return(MAE)
}
