#' @title reduceWiki
#' @description Returns all gene IDs of a single wikipathway of interest. This
#' function is recommended to be used after a signalling pathway of interest
#' is found.
#' @param MAE MultiAssayExperiment which will store the results of reduceWiki.
#' It is recommended to use the same MAE which stores the data from dloadGmt/
#' gmtEnsembl.
#' @param path_data Dataframe with wikipathway IDs, gene IDs and pathway names
#' from either the dloadGmt or gmtEnsembl functions. These will be found as
#' assays within the MAE used in the dloadGmt or gmtEnsembl functions.
#' @param stringWiki Full name of the wikipathway of interest. Make sure to
#' spell this correctly. Each wikipathway can only be added once to the same MAE
#' object.
#' @return A dataframe that only contains information about the wikipathway of
#' interest. Output will be stored as an assay in the input MAE.
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

    if (missing(MAE)) stop('MAE is missing. Add a MultiAssayExperiment to store results from reduceWiki. Please use the dloadGmt and/ or gmtEnsembl functions first.')

    if (missing(path_data)) stop('path_data is missing. Add dataframe with wikipathways IDs, gene IDs and wikipathway names. Please use dloadGmt and/ or gmtEnsembl first. Output of these functions will be found as assays within the MAE used in the dloadGmt and/ or gmtEnsembl functions.')

    if (missing(stringWiki)) stop('stringWiki is missing. Add name of chosen wikipathway e.g. "TGF-beta Signaling Pathway". Please make sure that a correctly spelled wikipathway name has been added.')

    # retrieve data from one pathway
    singlewiki <- path_data[which(path_data$name == stringWiki),]

    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = singlewiki)))

    # give new file unique name
    names(MAE2) <- stringWiki

    MAE <- suppressMessages(c(MAE, MAE2))
    return(MAE)
}
