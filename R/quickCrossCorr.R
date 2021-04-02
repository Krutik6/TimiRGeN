#' @title quickCrossCorr
#' @description Plots a cross-correlation plot to compare the miRNA and mRNA of
#' a selected pair. This is a useful test of the similarities between the the
#' two time series. It tracks movement of two time series relative to one
#' another to determine how well they match and at which point
#' the best match occurs.
#' @param filt_df Dataframe from the matrixFilter function.
#' @param pair Interger representing the pair to be explored.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data
#' @param scale TRUE or FALSE. Should data be scales. Default is FALSE. If
#' the correlation is based on Log2FC values scale should be TRUE.
#' @param Interpolation TRUE or FALSE. Should the whole time course be
#' interpolated over by a smooth spline? Default is FALSE. This is most useful
#' for longer and regular time courses.
#' @param timecourse If Iterpolation is TRUE, how many time points should be
#' interpolated over?
#' @return A cross correlation plot.
#' @usage quickCrossCorr(filt_df, pair, miRNA_exp, mRNA_exp, scale,
#' Interpolation, timecourse)
#' @export
#' @importFrom stats as.ts ccf
#' @importFrom scales rescale
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR[1:100,]
#'
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus')
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
#' quickCrossCorr(filt_df=MAE[[11]], pair=1, miRNA_exp=MAE[[9]],
#'               mRNA_exp=MAE[[10]],scale = FALSE, Interpolation = FALSE)
quickCrossCorr <- function(filt_df, pair, miRNA_exp, mRNA_exp, scale=FALSE,
                           Interpolation=FALSE, timecourse){

  Int <- pickPair(filt_df, pair, miRNA_exp, mRNA_exp, scale)


  if (Interpolation == TRUE) {

    if (missing(timecourse)) stop('timecourse is missing. How many time points to interpolate over? This should be the whole time course.')

    Int <- FreqProf::approxm(as.data.frame(Int), timecourse,

                             method = "spline")

  }else if (Interpolation == FALSE) {

    Int <- Int

  }

  DF <- as.data.frame(Int)

  M <- as.matrix(DF)

  M <- scales::rescale(M, to = c(1, 100))

  miR_tc <- as.numeric(M[,1])

  mRNA_tc <- as.numeric(M[,2])

  Y <- cbind(miR_tc, mRNA_tc)

  colnames(Y) <- colnames(DF[1:2])

  Z <- as.ts(Y)

  par(cex.main = 1.7)

  par(oma = c(1, 1, 1, 1))

  ccf(Z[,1], log(Z[,2]), ylab = "Cross-correlation", type = "correlation",

      main = paste0(names(DF)[1], ":", names(DF)[2],

                    " Cross-Correlation"), cex.lab=1.5, cex.axis=1.3, lwd =4)

}
