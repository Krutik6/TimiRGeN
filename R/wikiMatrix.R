#' @title wikiMatrix
#' @description Creates a matrix which shows how many genes from the data are
#'found in each wikipathway for a specific species.
#' @param MAE MultiAssayExperiment object.
#' @param ID_list list of gene names as strings.
#' @param wp_list list of list of wikipathways with gene names as strings.
#'
#' @return A matrix of timepoints and wikipathways.
#' @export
#' @usage wikiMatrix(MAE, ID_list, wp_list)
#'
#' @examples
#' MAE <- MultiAssayExperiment()
#' metadata(MAE)[["ID_list"]] <- e_list
#' metadata(MAE)[["w_list"]] <- w_list[1:10]
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
wikiMatrix <- function(MAE, ID_list, wp_list){

    if (missing(MAE)) stop('Add MAE object');
    if (missing(ID_list)) stop('Input list of list of genenames.');
    if (missing(wp_list)) stop('Input list of list of wikipathways.');

    wmat <- sapply(wp_list, function(x) {
        sapply(ID_list, function(y) sum(x %in% y))})

    L <- lapply(wp_list, function(x){length(x)})

    wmat2 <- as.data.frame(rbind(wmat, Total = unlist(L)))

    MAE2 <- suppressMessages(MultiAssayExperiment(list("wikimatrix" = wmat2)))
    MAE <- c(MAE, MAE2)

    return(MAE)
}
