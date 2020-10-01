#' @title matrixFilter
#' @description Filters out miR-mRNA interactions based on how many times an
#' interaction has been predicted and/ or validated. miR-mRNA interactions can
#' also be filtered by averaged correlations of DE values (log2fc or ave exp).
#' Negatively correlating miR-mRNA interactions can be filtered for, and degree
#' of correlation is also a filterable parameter.
#' @param MAE MultiAssayExperiment to store the output of matrixFilter. It is
#' recommended to use the same MAE which stores the results from
#' dataMiningMatrix.
#' @param miningMatrix A Large correlation matrix which has miR-mRNA
#' validation information from targetscans, mirdb and mirtarbase.
#' This is output from dataMiningMatrix, and should be stored as an assay within
#' the MAE used in the dataMiningMatrix function.
#' @param negativeOnly TRUE or FALSE. Should only negatively correlating
#' miR-mRNA interactions be retrieved? Default is TRUE.
#' @param predictedOnly TRUE or FALSE. Should only predicted interactions
#' should be retrieved? Default is TRUE.
#' @param threshold Integer from 0 to 3. How many databases should a miR-mRNA
#' interaction be found in? If predictedOnly = TRUE, then maximum threshold is
#' 2.
#' @param maxCor Number from -1 to 1. What is the highest average correlation
#' that is allowed? Default is 1. The lower the maxCor, the stricter the
#' filtering.
#' @return Filtered miR-mRNA interactions that are specific for a signalling
#' pathway of interest and the input data. Output will be stored as an assay
#' in the input MAE.
#' @export
#' @usage matrixFilter(MAE, miningMatrix, negativeOnly, predictedOnly,
#'                     threshold, maxCor)
#' @examples
#'Int_matrix <- data.frame(row.names = c("mmu-miR-320-3p:Acss1",
#'                                       "mmu-miR-27a-3p:Odc1"),
#'                         avecor = c(-0.9191653, 0.7826041),
#'                         miR = c("mmu-miR-320-3p", "mmu-miR-27a-3p"),
#'                         mRNA = c("Acss1", "Acss1"),
#'                         miR_Entrez = c(NA, NA),
#'                         mRNA_Entrez = c(68738, 18263),
#'                         TargetScan = c(1, 0),
#'                         miRDB = c(0, 0),
#'                         Predicted_Interactions = c(1, 0),
#'                         miRTarBase = c(0, 1),
#'                         Pred_Fun = c(1, 1))
#'
#' MAE <- MultiAssayExperiment()
#'
#' MAE <- matrixFilter(MAE, miningMatrix = Int_matrix, negativeOnly = TRUE,
#'                         threshold = 1, predictedOnly = FALSE)
matrixFilter <- function(MAE, miningMatrix, negativeOnly = TRUE,
                         predictedOnly = TRUE, threshold = 1, maxCor = 1){

    if (missing(MAE)) stop('
                           MAE is missing.
                           Add MAE. This will store the output of matrixFilter.
                           Please use dataMiningMatrix first.')

    if (missing(miningMatrix)) stop('
                                    miningMatrix is missing.
                                    Add large correlation matrix. Please
                                    use the dataMiningMatrix function first.
                                    The output of dataMiningMatrix should be
                                    found as an assay within the MAE used in
                                    the dataMiningMatrix function.')

    # miR-mRNA interactions with a - ave correlation
    if (negativeOnly == TRUE) {

        Filter1 <- miningMatrix[which(miningMatrix$avecor < 0),]

    } else Filter1 <- miningMatrix

    # miR-mRNA interactions from miRDB and TargetScans
    if (predictedOnly == TRUE) {

            Filter2 <- Filter1[which(
                            Filter1$Pred_only >= threshold),]

    } else Filter2 <- Filter1[which(Filter1$Pred_Fun >= threshold),]

    # Maximum average correlation
        Filter3 <- Filter2[which(Filter2$avecor <= maxCor),]

    MAE2 <- suppressMessages(MultiAssayExperiment(
                            list("Filtered_Int" = Filter3)))

    MAE <- c(MAE, MAE2)

return(MAE)
}
