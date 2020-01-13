#' @title ReduceWiki
#' @description Return gene names for a single wikipathway of interest.
#' @param path_data File with wikipathway - gene data.
#' @param stringWiki Full name of the wikipathway of interest.
#' @return A dataframe which only contains information about the wikipathway of
#'interest.
#' @export
#' @usage ReduceWiki(path_data, stringWiki = '')
#' @examples
#' path_data <- data.frame("wpid" = c(rep("WP571", 6)),
#' "gene" = c(16175, 12370,26419, 19249, 19645, 18479),
#' "name" = c(rep("Fas pathway and Stress induction of HSP regulation", 6)))
#' ReduceWiki(path_data = path_data,
#' stringWiki = 'Fas pathway and Stress induction of
#'HSP regulation') -> singlewiki
ReduceWiki <- function(path_data, stringWiki){
        if (missing(path_data)) stop('Input wikipathways - gene data.');
        if (missing(stringWiki)) stop('Input name of chosen wikipathway e.g.
        TGF-beta Signaling Pathway')
        path_data[which(path_data$name == stringWiki),] -> singlewiki
return(singlewiki)
}
