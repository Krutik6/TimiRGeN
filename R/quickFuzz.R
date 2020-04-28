#' @title quickFuzz
#' @description Plots Clusters in reference to the ClusterData dataframe. To be
#' used after CreateClusters and ClusterCheck. From here individual
#' wikipathways that show variance throughout a signalling network can be
#' taken for further miR-mRNA interaction analysis.
#' @param Mfuzzdata Output after CreateClusters. An Expressionset object.
#' @param Clusters Output after CreateClusters. A list.
#' @param W Should the plot be shown in a new window? Default is TRUE.
#' @return A plot of different clusters showing how the input data varies
#' across different wikipathways.
#' @export
#' @importFrom Mfuzz mfuzz.plot2
#' @usage quickFuzz(Mfuzzdata, Clusters, W)
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
#' MAE <- createClusters(MAE, method = "c",
#'                       percentMatrix = assay(MAE, 2),
#'                       noClusters = 2, variance = 0.99)
#'
#' quickFuzz(Mfuzzdata = experiments(MAE)[[4]],
#'           Clusters = metadata(MAE)[[3]], W = FALSE)
quickFuzz <- function(Mfuzzdata, Clusters, W = TRUE){

    lab <- colnames(Clusters$centers)
    Maxim <- max(unique(Clusters$cluster))
    a <- ceiling(Maxim/6)
    b <- ceiling(Maxim/a)
    mfuzz.plot2(eset = Mfuzzdata, cl = Clusters, mfrow = c(a, b),
                time.labels = lab, xlab = "Data points",
                ylab = "Normalised Common Genes between data and wikipathways",
                ax.col = "white", bg = "black", col.lab = "yellow",
                col.axis = "white", colo = "fancy", x11 = W)
}
