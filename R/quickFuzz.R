#' @title quickFuzz
#' @description Plots fuzzy clusters. Each different cluster created will
#' represent a different temporal behaviour. Depending on the data, more or
#' fewer cluster may be appropriate.
#' Use clusterCheck to influence this decision before moving onto quickFuzz.
#' Each line in a cluster represents a pathway. Pathways are divided by colour.
#' The more intense the colour of a line, the stronger they fit a particular
#' cluster / temporal behaviour.
#' Fuzzy clustering is a soft clustering approach where objects are not
#' divided into fixed clusters. Each pathway can exist in each cluster but
#' each pathway will differ on the degree to which they fit to each cluster.
#' Look into the clusterData dataframe created by createClusters to see this.
#' If a cluster peaks interest, continue to analysis of that cluster with the
#' returnCluster function.
#' @param Mfuzzdata A large ExpressionSet object which contain fuzzy clustering
#' data. This is output from the createClusters function. The Expressionset
#' object should be stored as an experiment in the MAE used in the
#' createClusters function.
#' @param Clusters A large list containing information about clusters,
#' statistics and phenodata. This is output from the createClusters function.
#' The list should be stored as metadata in the MAE used in the
#' createClusters function.
#' @param W TRUE or FALSE? Should the plot be shown in a new window? Default is
#' TRUE.
#' @param background Plot background colour. Default is black.
#' @param labelcol Plot labels colour. Default is yellow.
#' @param axiscol Plot axis labels colour. Default is white.
#' @param axisline Plot axis line colour. Default is white.
#' @param subcol Plot sub title colour. Default is yellow.
#' @param ylab y axis label. Default is "Genes found in data and pathway".
#' @return A plot of different clusters showing how the number of genes found
#' to be significant varies between the input data and wikipathways. These
#' variations are captured as temporal behaviours and are clustered.
#' @export
#' @importFrom Mfuzz mfuzz.plot2
#' @usage quickFuzz(Mfuzzdata, Clusters, W, background, labelcol, axiscol,
#' axisline, subcol, ylab)
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
#'                       percentMatrix = assay(MAE, 2),
#'                       noClusters = 2, variance = 0.99)
#'
#' quickFuzz(Mfuzzdata = experiments(MAE)[[4]],
#'           Clusters = metadata(MAE)[[3]], W = FALSE)
quickFuzz <- function(Mfuzzdata, Clusters, W = TRUE, background = "black",
                      labelcol = "yellow", axiscol = "white", axisline = "white",
                      subcol = 'yellow', ylab = "Genes found in data and pathway"){

    if (missing(Mfuzzdata)) stop('Mfuzzdata is missing. Add the ExpressionSet object created by createClusters. Please use createClusters first. Mfuzzdata should be stored as an experiment within the MAE used in the createClusters function.')

    if (missing(Clusters)) stop('Clusters is missing. Add the list created by createClusters. Please use createClusters first. Please use the createClusters function first. Clusters should be stored as metadata within the MAE used in the createClusters function.')

    lab <- colnames(Clusters$centers)

    Maxim <- max(unique(Clusters$cluster))

    a <- ceiling(Maxim/6)

    b <- ceiling(Maxim/a)

    Mfuzz::mfuzz.plot2(eset = Mfuzzdata,
                       cl = Clusters,
                       mfrow = c(a, b),
                       time.labels = lab,
                       xlab = "Time points",
                       ylab = ylab,
                       ax.col = axisline,
                       bg = background,
                       col.lab = labelcol,
                       col.axis = axiscol,
                       colo = "fancy",
                       col.main = subcol,
                       x11 = W)
}
