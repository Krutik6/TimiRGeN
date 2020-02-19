#' @title ClusterCheck
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
#' @usage ClusterCheck(Clusters, W)
#' @examples
#' library(Mfuzz)
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'e_list -> MAE@metadata$elist
#'w_list[1:10] -> MAE@metadata$wlist
#'WikiMatrix(e_list = MAE@metadata$elist, 
#'wp_list = MAE@metadata$wlist) -> MAE@ExperimentList$Wmat
#'
#'TurnPercent(wikiMatrix = MAE@ExperimentList$Wmat,
#'rowInt = 4) -> MAE@ExperimentList$Pmat
#'
#'CreateClusters(method = "c", MAE, Percent_matrix = MAE@ExperimentList$Pmat,
#'no.clusters = 2, Variance = 0.99) -> MAE
#' 
#' ClusterCheck(Clusters = MAE@metadata$Clusters, W = FALSE)
ClusterCheck <- function(Clusters, W = TRUE){
dev.new <- NULL
O <- overlap(cl = Clusters)
if (W == TRUE) {
dev.new()
overlap.plot(Clusters, overlap=O, thres=0.05)
} else overlap.plot(Clusters, overlap=O, thres=0.05)
}
