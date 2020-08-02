#' @title wikiMatrix
#' @description Creates a matrix which shows how many genes from the input mRNA
#' and miRNA data are found in each wikipathway for a specific species. The
#' matrix is stored as an assay in a MAE.
#' @param MAE MultiAssayExperiment which will have the results of
#' wikiMatrix stored  within it as an assay. It is recommended that the
#' MAE object which stores output from the wikiList to be used here.
#' @param ID_list List of entrez gene IDs or ensembl gene IDs stored as strings.
#' This is the output of the eNames function and is stored as metadata in the
#' MAE used in the eNames function.
#' @param wp_list List of list containing wikipathways with entrez gene IDs or
#' ensembl IDs as strings. This is the output of the wikiList function and can
#' be found as metadata in the MAE used in the wikiList function.
#' @return A matrix of gene found in each time points for each wikipathways.
#' @export
#' @usage wikiMatrix(MAE, ID_list, wp_list)
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["ID_list"]] <- e_list_mouse
#'
#' metadata(MAE)[["w_list"]] <- w_list_mouse[1:10]
#'
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
wikiMatrix <- function(MAE, ID_list, wp_list){

    if (missing(MAE)) stop('Add MAE object which will have the results of
                           wikiMatrix added. Please use the eNames and
                           wikiMatrix functions first.')

    if (missing(ID_list)) stop('Add list of list of gene IDs  Please
                               use the eNames function before using wikiMatrix.
                               Output of eNames should be stored as metadata
                               within the MAE used in the eNames function.')

    if (missing(wp_list)) stop('Add list of list of wikipathways
                               and associated entrezgene or ensembl IDs.
                               Please use the wikiList function before using
                               wikiMatrix. Output of wikiList should be stored
                               as metadata within the MAE used in the
                               wikiList function.')

    # Convert list into matrix
    wmat <- sapply(wp_list, function(x) {
        sapply(ID_list, function(y) sum(x %in% y))})

    L <- lapply(wp_list, function(x){length(x)})

    wmat2 <- as.data.frame(rbind(wmat, Total = unlist(L)))

    MAE2 <- suppressMessages(MultiAssayExperiment(list("wikimatrix" = wmat2)))

    MAE <- c(MAE, MAE2)

    return(MAE)
}
