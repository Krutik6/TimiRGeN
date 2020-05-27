#' @title returnCluster
#' @description Will retrieve information about which wikipathways fitted best
#' with a specific cluster.
#' @param MAE MultiAssayExperiment object.
#' @param clusterData Output after CreateClusters. A list.
#' @param whichCluster Integer should correspond to the order of clusters
#' displayed.
#' @param fitCluster How well should the wikipathways fit into the this
#' cluster? Integer from 0-1. Default is 0.99.
#' @return A dataframe of which pathways corresponded best with the chosen
#' dynamics seen in the selected cluster.
#' @export
#' @usage returnCluster(MAE, clusterData, whichCluster, fitCluster)
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list
#'
#' metadata(MAE)[["w_list"]] <- w_list[1:10]
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

    if (missing(MAE)) stop('Add MAE object.')

    if (missing(clusterData)) stop('Add clusterdata data frame from
                                   createClusters function. Will be in
                                   Experimentlist.')

    if (missing(whichCluster)) stop('Add interger which represent which
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
