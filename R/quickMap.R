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
#' library(org.Mm.eg.db)
#' miR <- mm_miR[1:100,]
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus')
#'
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 1), dataType = 'Log2FC',
#'                      genes_ID = assay(MAE, 3),
#'                       idColumn = 'GENENAME',
#'                       name = "miRNA_log2fc")
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                      genes_ID = assay(MAE, 7),
#'                      idColumn = 'GENENAME',
#'                      name = "mRNA_log2fc")
#'
#' Filt_df <- data.frame(row.names = c("mmu-miR-145a-3p:Adamts15",
#'                                    "mmu-miR-146a-5p:Acy1"),
#'                      corr = c(-0.9191653, 0.7826041),
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
#' MAE <- matrixFilter(MAE, miningMatrix = Filt_df, negativeOnly = FALSE,
#'                    threshold = 1, predictedOnly = FALSE)
#'
#' quickMap(filt_df = MAE[[11]], numpairs = 2)
quickMap <- function(filt_df, numpairs){

  if (missing(filt_df)) stop('filt_df is missing. Add the assay/ dataframe which is the output from matrixFiler.')

  if (missing(numpairs)) stop('numpairs is missing. Add an integer greater than 1.')

  Ranks <- filt_df[,c(1,2,3)][order(filt_df$corr,

                                    decreasing = FALSE),]

  Ranks$miR <- Ranks$mRNA <- NULL

  Ranks <- as.matrix(Ranks)

  Ranks <- cbind(Ranks[,1], Ranks[,1])

  par(oma = c(0, 2, 0, 7))

  par(cex.main=1.5)

  Ranks[order(Ranks)]

  heatmap.2(Ranks[1:numpairs,],

            trace = "n",

            labCol = "",

            ylab = "",

            Colv = FALSE,

            dendrogram = "none",

            labRow = rownames(Ranks[1:numpairs,]),

            main = "miR-mRNA pairs",

            cexRow = 1.15,

            key.title = "",

            Rowv=FALSE,

            rowsep=c(1:numpairs),

            sepcolor="white",

            key.xlab = "Corr",

            key.ylab = "")

}
