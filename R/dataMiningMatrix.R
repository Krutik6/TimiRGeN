#' @title dataMiningMatrix
#' @description Mines out predicted/ functional interactions which correspond
#' between miR-mRNA interactions found in Targetscans, miRDB, miRTarBase and
#' the interactions from the miR-mRNA correlation matrix.
#' If a database cannot be downloaded, dataMiningMatrix can be used regardless,
#' but it is recommended to download all three databases.
#' @param MAE MultiAssayExperiment which will store the output of
#' dataMiningMatrix. It is recommended to use the MAE object which stores output
#' from the mirMrnaInt function.
#' @param corrTable Correlation matrix of interactions between the mRNAs from
#' a pathway of interest and miRNA data. This is created by the mirMrnaInt
#' function and should be stored as an assay within the MAE used in the
#' mirMrnaInt function.
#' @param targetscan Species specific miR-mRNA interactions predicted by
#' targetscans. This is the output from the dloadTargetscan function.
#' It should be stored as an assay within the MAE used in the dloadTargetscan
#' function. If this data cannot be downloaded, dataMiningMatrix can be run
#' without it.
#' @param mirdb Species specific miR-mRNA interactions predicted by miRDB.
#' This is the output from the dloadMirdb function. It should
#' stored as an assay within the MAE used in the dloadMirdb function. If
#' this data cannot be downloaded, dataMiningMatrix can be run without it.
#' @param mirtarbase Species specific miR-mRNA interactions which are
#' functionally curated by mirtarbase. This is the output from the
#' dloadMirtarbase function. It should be stored as an assay within the MAE used
#' in the dloadMirtarbase function. If this data cannot be downloaded,
#' dataMiningMatrix can be run without it.
#' @return A matrix which cross references the occurrences of miR-mRNA
#' interactions between databases and the given data. Output will be stored as
#' an assay in the input MAE.
#' @export
#' @usage dataMiningMatrix(MAE, corrTable, targetscan , mirdb, mirtarbase)
dataMiningMatrix <- function(MAE, corrTable, targetscan , mirdb,
                             mirtarbase){

  if (missing(MAE)) stop('MAE is missing. Add MAE. This will store the output of dataMiningMatrix. Please use the mirMrnaInt, dloadTargetscan, dloadMirdb and dloadMirtarbase functions first.')

  if (missing(corrTable)) stop('corrTable is missing. Add matrix of miR-mRNA interactions and correlations. Please use the mirMrnaInt function first. The output of mirMrnaInt should be found as an assay within the MAE used in the mirMrnaInt function.')

  TargetScan <- miRDB <- Predicted_Interactions <- miRTarBase <- NULL

  X <- corrTable

  X$Pred_Fun <- X$miRTarBase <- X$Pred_only <- X$miRDB <-
                X$TargetScan <-numeric(nrow(X))

  if(!missing(targetscan)) X$TargetScan <- as.integer(
                                              rownames(X) %in% targetscan[[1]])

  if(!missing(mirdb)) X$miRDB <- as.integer(rownames(X) %in% mirdb[[1]])

  X$Pred_only <- X$TargetScan + X$miRDB

  if(!missing(mirtarbase)) X$miRTarBase <- as.integer(
                                              rownames(X) %in% mirtarbase[[1]])

  X$Pred_Fun <- X$Pred_only + X$miRTarBase

  # Add to MAE object
  MAE2 <- suppressMessages(MultiAssayExperiment(list("Int_Matrix" = X)))

  MAE <- c(MAE, MAE2)

  return(MAE)
}
