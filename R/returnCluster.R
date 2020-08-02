#' @title returnCluster
#' @description Will retrieve information about which wikipathways fitted best
#' with a specific cluster and store this as an assay within a MAE.
#' @param MAE MultiAssayExperiment. Results from returnCluster will
#' be stored within this MAE. It is recommended to use the same MAE which
#' stores output from createClusters.
#' @param clusterData A dataframe which contains cluster-pathway fit scores
#' and is stored as an assay within the MAE used in the createClusters function.
#' @param whichCluster Integer should corresponds to the cluster of interest.
#' @param fitCluster How well should the pathways fit into the this
#' cluster? Integer from 0-1. Default is 0.99.
#' @return A dataframe of which pathways corresponded best with the chosen
#' dynamics seen in the selected cluster.
#' @export
#' @usage returnCluster(MAE, clusterData, whichCluster, fitCluster)
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
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
#'
#' MAE <- createClusters(MAE, method = "c",
#'                     percentMatrix = assay(MAE, 2),
#'                     noClusters = 2, variance = 0.99)
#'
#' MAE <- returnCluster(MAE, clusterData = assay(MAE, 3), whichCluster = 1,
#'                      fitCluster = 0.5)
returnCluster <-function(MAE, clusterData, whichCluster, fitCluster = 0.99){

    if (missing(MAE)) stop('Add MAE object. Results from returnCluster
                           will be stored within given MAE. Please
                           use createClusters first.')

    if (missing(clusterData)) stop('Add dataframe which has cluster-pathway
                                   fit scores. Please use the createClusters
                                   function. clusterData should be
                                   stored as an assay within the MAE used in
                                   the createClusters function.')

    if (missing(whichCluster)) stop('Add interger which represents which
                                    cluster is of interest.')

    X <- as.data.frame(clusterData)


    # Remove pathways from cluster which have a fit score lower than fitCluster
    singlecluster <- X[which(X[,whichCluster] > fitCluster),]

    # Make a unique name for the files
    a <- "Cluster:"

    cl <- whichCluster

    b <- "_fit:"

    fit <- fitCluster

    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = singlecluster)))

    # Paste unique name
    names(MAE2) <- paste0(a, cl, b, fit)

    MAE <- c(MAE, MAE2)
return(MAE)
}
