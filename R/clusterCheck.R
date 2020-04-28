#' @title clusterCheck
#' @description Put data from CreateClusters through this funciton to check
#' if an appropriate number of clusters have been created. The closer the
#' circles are to one another the more likely that they should belong
#' to the same cluster. Read more about Mfuzz for a greater understanding
#' https://bioconductor.org/packages/release/bioc/html/rWikiPathways.html.
#' @param Clusters Outout from CreateClusters function. A list of clusters,
#' statistics and phenodata.
#' @param W Should the plot be shown in a new window? Default is TRUE.
#' @return A PCAplot showing distance of clusters.
#' @export
#' @importFrom grDevices dev.new
#' @importFrom Mfuzz overlap overlap.plot
#' @usage clusterCheck(Clusters, W)
#' @examples
#' MAE <- MultiAssayExperiment()
#' metadata(MAE)[["e_list"]] <- e_list
#' metadata(MAE)[["w_list"]] <- w_list[1:10]
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
#'
#' MAE <- turnPercent(MAE = MAE,
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
#'
#'MAE <- createClusters(MAE, method = "c",
#'                     percentMatrix = assay(MAE, 2),
#'                     noClusters = 2, variance = 0.99)
#'
#' clusterCheck(Clusters = metadata(MAE)[[3]], W = FALSE)
clusterCheck <- function(Clusters, W = TRUE){
    if (missing(Clusters)) stop('Add list of Clusters from CreateClusters');

    dev.new <- NULL
    #Check Clusters using a PCA plot
    O <- overlap(cl = Clusters)
    if (W == TRUE) {
        dev.new()
        overlap.plot(Clusters, overlap=O, thres=0.05)
    } else overlap.plot(Clusters, overlap=O, thres=0.05)
}
