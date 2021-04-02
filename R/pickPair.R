#' @title pickPair
#' @description Internal function to be used for quickTC and quickCrossCorr
#' @param filt_df Dataframe from the matrixFilter function.
#' @param pair Interger representing the pair to be explored.
#' @param miRNA_exp miRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param mRNA_exp mRNA data from using the diffExpressRes function on miRNA
#' data.
#' @param scale TRUE od FALSE. Sould the data be scaled? Default is FALSE.
#' @return A time course of a miRNA:mRNA pair of interest, plotted along time.
#' @noRd
#' @usage pickPair(filt_df, pair, miRNA_exp, mRNA_exp, scale)
pickPair <- function(filt_df, pair, miRNA_exp, mRNA_exp, scale = FALSE){

  if (missing(filt_df)) stop('filt_df is missing. Add assay/ dataframe created by the matrixFilter function.')

  if (missing(miRNA_exp)) stop('miRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on miRNA expression data/ DE data.')

  if (missing(mRNA_exp)) stop('mRNA_exp is missing. Add assay/ dataframe created by the diffExpressRes function used on mRNA expression data/ DE data.')

  Ranks <- filt_df[,c(1,2,3)][order(filt_df$corr, decreasing = FALSE),]

  miRNA_target <- Ranks[[2]][pair]

  mRNA_target <- Ranks[[3]][pair]

  miR_tc <- miRNA_exp[which(rownames(miRNA_exp) == miRNA_target),]

  mRNA_tc <- mRNA_exp[which(rownames(mRNA_exp) == mRNA_target),]

  miR_tc$ID <- NULL

  mRNA_tc$ID <- NULL

  times <- as.integer(gsub(colnames(miR_tc), pattern = "[^0-9.-]",

                           replacement = ""))

  times <- as.numeric(times)

  miR_tc <- as.numeric(miR_tc)

  mRNA_tc <- as.numeric(mRNA_tc)

  if (scale == TRUE) {

    miR_tc <- scale(miR_tc)

    mRNA_tc <- scale(mRNA_tc)

    miR_tc <- as.numeric(miR_tc)

    mRNA_tc <- as.numeric(mRNA_tc)

  } else if(scale == FALSE) {

    miR_tc <- miR_tc

    mRNA_tc <- mRNA_tc

  }
  Y <- cbind(miR_tc, mRNA_tc, times)

  colnames(Y) <- c(miRNA_target, mRNA_target, "Time")

  return(Y)
}
