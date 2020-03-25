#' @title Quickfuzz
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
#' @usage Quickfuzz(Mfuzzdata, Clusters, W)
#' @examples
#' library(Mfuzz)
#' library(MultiAssayExperiment)
#' MAE <- Clusters
#' Quickfuzz(Mfuzzdata = experiments(MAE)[[2]],
#'           Clusters = metadata(MAE)[[1]], W = FALSE)
Quickfuzz <- function(Mfuzzdata, Clusters, W = TRUE){
    
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
