#' @title reduceWiki
#' @description Return gene names for a single wikipathway of interest and
#' store this as an assay within a MAE.
#' @param MAE MultiAssayExperiment object which will have the results of
#' reduceWiki stored within as an assay. It is recommended to use the
#' same MAE object which stores the data from dloadGmt/ gmtEnsembl. Each
#' wikipathway can only be added once to the same MAE object.
#' @param path_data Dataframe with wikipathway IDs, gene IDs and pathway names
#' from dloadGmt or gmtEnsembl. This will be found as an assay in the MAE
#' used in the dloadGmt or gmtEnsembl functions.
#' @param stringWiki Full name of the wikipathway of interest. Make sure to
#' spell it correctly.
#' @return A dataframe which only contains information about the wikipathway of
#'interest.
#' @export
#' @usage reduceWiki(MAE, path_data, stringWiki = '')
#' @examples
#' MAE <- MultiAssayExperiment(list(path_data = data.frame(
#'                         "wpid" = c(rep("WP571", 6)),
#'                         "gene" = c(16175, 12370,26419, 19249, 19645, 18479),
#'                         "name" = c(rep("Fas pathway and Stress induction of HSP regulation",
#'                         6)))))
#'
#' MAE <- reduceWiki(MAE, path_data = assay(MAE, 1),
#'                  stringWiki = 'Fas pathway and Stress induction of HSP regulation')
reduceWiki <- function(MAE, path_data, stringWiki){

    if (missing(MAE)) stop('Add a MultiAssayExperiment object to store results
                           from reduceWiki. Have you used dloadGmt and/ or
                           gmtEnsembl?')

    if (missing(path_data)) stop('Input dataframe with wikipathways IDs,
                                 gene IDs and wikipathway names from.
                                 Please use dloadGmt and/ or gmtEnsembl first.
                                 Output of these functions will be found as
                                 assays within the MAE used in the dloadGmt
                                 and/ or gmtEnsembl functions.')

    if (missing(stringWiki)) stop('Input name of chosen wikipathway e.g.
                                  TGF-beta Signaling Pathway. Please make sure
                                  you have added a correct wikipathway name.')

    # retrieve data from one pathway
    singlewiki <- path_data[which(path_data$name == stringWiki),]

    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = singlewiki)))

    # give new file unique name
    names(MAE2) <- stringWiki

    MAE <- suppressMessages(c(MAE, MAE2))
    return(MAE)
}
