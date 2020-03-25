#' @title MatrixFilter
#' @description Filters out miR-mRNA interactions based on how many times an
#'interaction has been predicted and or validated.
#' @param MAE MultiAssayExperiment.
#' @param miningMatrix Output from DataMiningMatrix.
#' @param negativeOnly TRUE or FALSE. Should only negatively correlating
#' miR-mRNA interactions be given? Default is TRUE.
#' @param predictedOnly TRUE or FALSE. Should only predicted interactions be
#' displayed or predicted and functionally tested? Default is TRUE.
#' @param threshold Integer from 0 to 3. If Predicted only max should be 2.
#' @param maxCor What is the highest average correlation to be filtered for.
#' Default is 1.
#' @return A filtered amount of mRNA-miR interactions which are found in X
#' number of databases. This will be specific for the input data and the
#' wikipathway of choice.
#' @export
#' @usage MatrixFilter(MAE, miningMatrix, negativeOnly, predictedOnly, 
#'                     threshold, maxCor)
#' @examples
#'library(MultiAssayExperiment)
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
#' MAE <- MultiAssayExperiment()
#' MAE <- MatrixFilter(MAE, miningMatrix = Int_matrix, negativeOnly = TRUE,
#'                         threshold = 3, predictedOnly = FALSE)
MatrixFilter <- function(MAE, miningMatrix, negativeOnly = TRUE,
                         predictedOnly = TRUE, threshold = 1, maxCor = 1){
    
    if (missing(MAE)) stop('Add MAE object');
    if (missing(miningMatrix)) stop('Input miningmatrix from
                                     DataMiningMatrix function');
    
    if (negativeOnly == TRUE) {
        Filter1 <- miningMatrix[which(miningMatrix$avecor < 0),]
        
    } else Filter1 <- miningMatrix 
    
    if (predictedOnly == TRUE) {
            Filter2 <- Filter1[which(
                            Filter1$Predicted_Interactions >= threshold),]
            
    } else Filter2 <- Filter1[which(Filter1$Pred_Fun >= threshold),]
        Filter3 <- Filter2[which(Filter2$avecor <= maxCor),]
        
    MAE2 <- suppressMessages(MultiAssayExperiment(
                            list("Filtered_Int" = Filter3)))
    MAE <- c(MAE, MAE2)
    
return(MAE)
}
