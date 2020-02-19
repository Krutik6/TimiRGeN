#' @title ReturnCluster
#' @description Will retrive information about which wikipathways fitted best
#' with a specific cluster.
#' @param ClusterData Output after CreateClusters. A list.
#' @param which.cluster Integer should correspond to the order of clusters
#' displayed.
#' @param fit.cluster How well should the wikipathways fit into the this
#' cluster? Integer from 0-1. Default is 0.99.
#'
#' @return A dataframe of which pathways corresponded best with the chosen
#' dynamics seen in the selected cluster.
#' @export
#'
#' @usage ReturnCluster(ClusterData, which.cluster, fit.cluster)
#' @examples
#' library(Mfuzz)
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'e_list -> MAE@metadata$elist
#'w_list[1:10] -> MAE@metadata$wlist
#'WikiMatrix(e_list = MAE@metadata$elist, 
#' wp_list = MAE@metadata$wlist) -> MAE@ExperimentList$Wmat
#'
#'TurnPercent(wikiMatrix = MAE@ExperimentList$Wmat,
#'rowInt = 4) -> MAE@ExperimentList$Pmat
#'
#'CreateClusters(method = "c", MAE =  MAE,
#'Percent_matrix = MAE@ExperimentList$Pmat,
#'no.clusters = 2, Variance = 0.99) -> MAE
#' 
#' ReturnCluster(ClusterData = MAE@ExperimentList$ClusterData, 
#' which.cluster = 1, fit.cluster = 0.9) -> MAE@ExperimentList$CLUST1
ReturnCluster <-function(ClusterData, which.cluster, fit.cluster = 0.99){
X <- as.data.frame(ClusterData)
singlecluster <- X[which(X[,which.cluster] > fit.cluster),]
return(as.data.frame(singlecluster))
}
