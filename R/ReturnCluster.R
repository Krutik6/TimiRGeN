#' @title ReturnCluster
#' @description Will retrive information about which wikipathways fitted best
#' with a specific cluster.
#' @param MAE MultiAssayExperiment object.
#' @param clusterData Output after CreateClusters. A list.
#' @param whichCluster Integer should correspond to the order of clusters
#' displayed.
#' @param fitCluster How well should the wikipathways fit into the this
#' cluster? Integer from 0-1. Default is 0.99.
#'
#' @return A dataframe of which pathways corresponded best with the chosen
#' dynamics seen in the selected cluster.
#' @export
#'
#' @usage ReturnCluster(MAE, clusterData, whichCluster, fitCluster)
#' @examples
#' library(Mfuzz)
#' library(MultiAssayExperiment)
#' MAE <- MultiAssayExperiment()
#' metadata(MAE)[["e_list"]] <- e_list
#' metadata(MAE)[["w_list"]] <- w_list[1:10]
#' MAE <- WikiMatrix(MAE, ID_list = metadata(MAE)[[1]], 
#'                   wp_list = metadata(MAE)[[2]])
#'  
#' MAE <- TurnPercent(MAE = MAE, 
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
#' 
#'MAE <- CreateClusters(MAE, method = "c", 
#'                    percentMatrix = assay(MAE, 2),
#'                    noClusters = 2, variance = 0.99)
#' 
#' MAE <- ReturnCluster(MAE, clusterData = assay(MAE, 3), whichCluster = 1, 
#'                      fitCluster = 0.5)
ReturnCluster <-function(MAE, clusterData, whichCluster, fitCluster = 0.99){
    X <- as.data.frame(clusterData)
    singlecluster <- X[which(X[,whichCluster] > fitCluster),]
    # Make a unique name for the files
    a <- "Cluster:" 
    cl <- whichCluster
    b <- "_fit:"
    fit <- fitCluster
    
    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = singlecluster)))
    names(MAE2) <- paste0(a, cl, b, fit)
    
    MAE <- c(MAE, MAE2)
return(MAE)
}
