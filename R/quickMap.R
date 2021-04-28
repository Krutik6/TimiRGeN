#' @title quickMap
#' @description Generates a heatmap of all miRNA:mRNA binding pairs that
#' have been filtered. Pairs are ranked by decreasing correlation.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param numpairs Number of pairs to plot. Must be an integer more than 1.
#' @return Heatmap of miRNA-mRNA pairs.
#' @export
#' @usage quickMap(filt_df, numpairs)
#' @importFrom gplots heatmap.2
#' @examples
#' Filt_df <- data.frame(row.names = c("mmu-miR-145a-3p:Adamts15",
#'                                    "mmu-miR-146a-5p:Acy1"),
#'                      corr = c(-0.9191653, -0.7826041),
#'                      miR = c("mmu-miR-145a-3p", "mmu-miR-146a-5p"),
#'                      mRNA = c("Adamts15", "Acy1"),
#'                      miR_Entrez = c(387163, NA),
#'                      mRNA_Entrez = c(235130, 109652),
#'                      TargetScan = c(1, 0),
#'                      miRDB = c(0, 0),
#'                      Predicted_Interactions = c(1, 0),
#'                      miRTarBase = c(0, 1),
#'                      Pred_Fun = c(1, 1))
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- matrixFilter(MAE, miningMatrix = Filt_df, negativeOnly = FALSE,
#'                    threshold = 1, predictedOnly = FALSE)
#'
#' quickMap(filt_df = MAE[[1]], numpairs = 2)
quickMap <- function(filt_df, numpairs){

  if (missing(filt_df)) stop('filt_df is missing. Add the assay/ dataframe which is the output from matrixFiler.')

  if (missing(numpairs)) stop('numpairs is missing. Add an integer greater than 1.')

  tracecol <- NULL

  Ranks <- filt_df[,c(1,2,3)][order(filt_df$corr, decreasing = FALSE),]

  Ranks$miR <- Ranks$mRNA <- NULL

  Ranks <- as.matrix(Ranks)

  Ranks <- cbind(Ranks[,1], Ranks[,1])

  par(oma = c(0, 2, 0, 12))

  par(cex.main=1.5)

  Ranks[order(Ranks)]

  my_palette <- grDevices::colorRampPalette(c("red", "yellow", "white", "cyan",
                                              "blue"))(n = 200)

  x <- seq(1:numpairs)

  heatmap.2(Ranks[1:numpairs,],
            trace = "n",
            labCol = "",
            ylab = "",
            Colv = FALSE,
            dendrogram = "none",
            labRow = paste0(x, " : ", rownames(Ranks[1:numpairs,])),
            main = "miR-mRNA pairs",
            cexRow = 1.2,
            key.title = "",
            Rowv=FALSE,
            rowsep=c(1:numpairs),
            sepcolor="white",
            key.xlab = "Corr",
            key.ylab = "",
            col = my_palette,
            density.info=c("none"),
            denscol=tracecol)

}
