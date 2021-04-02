#' @title clusterCheck
#' @description clusterCheck creates a PCA plot using functions from Mfuzz. This
#' will help to indicate if an appropriate number of clusters have been created.
#' The closer the circles are to one another the more likely that they should
#' belong to the same cluster. Read more about Mfuzz here
#' https://bioconductor.org/packages/release/bioc/html/rWikiPathways.html.
#' @param Clusters A large list of clusters, statistics and phenodata. This
#' will be stored as metadata within the MAE used in the createClusters
#' function.
#' @param W TRUE or FALSE. Should the plot be shown in a new window? Default is
#' FALSE.
#' @return A PCAplot showing distance of clusters.
#' @export
#' @importFrom grDevices dev.new
#' @importFrom Mfuzz overlap overlap.plot
#' @usage clusterCheck(Clusters, W)
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
#'                    wikiMatrix = assay(MAE, 1))
#'
#' MAE <- createClusters(MAE, method = "c",
#'                     percentMatrix = assay(MAE, 2),
#'                     noClusters = 2, variance = 0.99)
#'
#' clusterCheck(Clusters = metadata(MAE)[[3]], W = FALSE)
clusterCheck <- function(Clusters, W = FALSE){

    if (missing(Clusters)) stop('Clusters is missing. Please use the createClusters or createClusters2 function first. The Clusters list created from this will be stored as metadata within the MAE used in the createClusters function.')

    dev.new <- NULL

    #Check Clusters using a PCA plot
    O <- Mfuzz::overlap(cl = Clusters)

    if (W == TRUE) {

        grDevices::dev.new()

        Mfuzz::overlap.plot(Clusters, overlap=O, thres=0.05)

    } else overlap.plot(Clusters, overlap=O, thres=0.05)
}
