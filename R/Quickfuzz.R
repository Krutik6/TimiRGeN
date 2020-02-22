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
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'e_list -> MAE@metadata$elist
#'w_list[1:10] -> MAE@metadata$wlist
#'WikiMatrix(e_list = MAE@metadata$elist, 
#' wp_list = MAE@metadata$wlist) -> MAE@ExperimentList$Wmat
#'
#'TurnPercent(wikiMatrix = MAE@ExperimentList$Wmat,
#'rowInt = 4) -> MAE@metadata$Pmat
#'
#'CreateClusters(method = "c", MAE =  MAE,
#'Percent_matrix = MAE@metadata$Pmat,
#' no.clusters = 2, Variance = 0.99) -> MAE
#' 
#' Quickfuzz(Mfuzzdata = MAE@ExperimentList$MfuzzData,
#' Clusters = MAE@metadata$Clusters, W = FALSE)
Quickfuzz <- function(Mfuzzdata, Clusters, W = TRUE){
lab <- colnames(Clusters$centers)
Maxim <- max(unique(Clusters$cluster))
ceiling(Maxim/6) -> a
ceiling(Maxim/a) -> b
mfuzz.plot2(eset = Mfuzzdata, cl = Clusters, mfrow = c(a, b),
time.labels = lab, xlab = "Data points",
ylab = "Normalised Common Genes between data and wikipathways",
ax.col = "white",
bg = "black", col.lab = "yellow", col.axis = "white",
colo = "fancy", x11 = W)
}
