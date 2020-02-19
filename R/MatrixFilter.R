#' @title MatrixFilter
#' @description Filters out miR-mRNA interactions based on how many times an
#'interaction has been predicted and or validated.
#' @param miningMatrix Output from DataMiningMatrix.
#' @param NegativeOnly TRUE or FALSE. Should only negatively correlating
#' miR-mRNA interactions be given? Default is TRUE.
#' @param PredictedOnly TRUE or FALSE. Should only predicted interactions be
#' displayed or predicted and functionally tested? Default is TRUE.
#' @param THRESHOLD Integer from 0 to 3. If Predicted only max should be 2.
#' @param maxCor What is the highest average correlation to be filtered for.
#' Default is 1.
#' @return A filtered amount of mRNA-miR interactions which are found in X
#' number of databases. This will be specific for the input data and the
#' wikipathway of choice.
#' @export
#' @usage MatrixFilter(miningMatrix, NegativeOnly, PredictedOnly, THRESHOLD,
#' maxCor)
#' @examples
#'Int_matrix <- data.frame(row.names = c("mmu-miR-320-3p:Acss1",
#' "mmu-miR-27a-3p:Odc1"),
#' avecor = c(-0.9191653, 0.7826041),
#' miR = c("mmu-miR-320-3p", "mmu-miR-27a-3p"),
#' mRNA = c("Acss1", "Acss1"),
#' miR_Entrez = c(NA, NA),
#' mRNA_Entrez = c(68738, 18263),
#' TargetScan = c(1, 0),
#' miRDB = c(0, 0),
#' Predicted_Interactions = c(1, 0),
#' miRTarBase = c(0, 1),
#' Pred_Fun = c(1, 1))
#' MatrixFilter(miningMatrix = Int_matrix ) -> Filt_df
#' MatrixFilter(miningMatrix = Int_matrix, NegativeOnly = TRUE,
#' THRESHOLD = 3, PredictedOnly = FALSE) -> Filt_df
MatrixFilter <- function(miningMatrix, NegativeOnly = TRUE,
PredictedOnly = TRUE, THRESHOLD = 1, maxCor = 1){
if (missing(miningMatrix)) stop('Input miningmatrix from
DataMiningMatrix function');
if (NegativeOnly == TRUE) {
miningMatrix[which(miningMatrix$avecor < 0),] -> Filter1
} else miningMatrix -> Filter1
if (PredictedOnly == TRUE) {
Filter1[which(Filter1$Predicted_Interactions >= THRESHOLD),] -> Filter2
} else Filter1[which(Filter1$Pred_Fun >= THRESHOLD),] -> Filter2
Filter2[which(Filter2$avecor <= maxCor),] -> Filter3
return(Filter3)
}
