#' @title dataMiningMatrix
#' @description Mines out predicted/ functional interactions which correspond
#'between Targetscans, miRDB, miRTarBase and the interactions data frame.
#' @param MAE MultiAssayExperiment.
#' @param corrTable Correlation matrix between the mRNAs from the
#' wikipathway of interest and miRNA data.
#' @param targetscan Targetscans predicted interactions.
#' @param mirdb miRDB predcited interactions.
#' @param mirtarbase miRTarBase functional interactions.
#' @return A matrix of information which cross references the occurance of
#'microRNA-mRNA interactions between databases and the given data.
#' @export
#' @importFrom stringr %>%
#' @importFrom dplyr mutate
#' @usage dataMiningMatrix(MAE, corrTable, targetscan , mirdb, mirtarbase)
dataMiningMatrix <- function(MAE, corrTable, targetscan, mirdb,
                             mirtarbase){

    if (missing(corrTable)) stop('Input interactions
                                 (miR-mRNA(wikipathway)) dataframe.');
    if (missing(targetscan)) stop('Input filtered targetscans file,
                                   interactions column.');
    if (missing(mirdb)) stop('Input filtered targetscans file,
                              interactions column.');
    if (missing(mirtarbase)) stop('Input filtered targetscans file,
                                    interactions column.');

    TargetScan <- miRDB <- Predicted_Interactions <- miRTarBase <- NULL
    # Are miR-mRNA interactions found in Targetscans?
    df2 <- corrTable %>% mutate(TargetScan = as.integer(
                                      rownames(corrTable) %in% targetscan[[1]]))
    rownames(df2) <- rownames(corrTable)
    # Are miR-mRNA interactions found in miRDB?
    df3 <- df2 %>% mutate(miRDB = as.integer(rownames(df2) %in% mirdb[[1]]))
    rownames(df3) <- rownames(corrTable)
    # Which occur in both?
    df4 <- df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB)
    rownames(df4) <- rownames(corrTable)
    # Are miR-mRNA interactions found in mirtarbase?
    df5 <- df4 %>% mutate(miRTarBase = as.integer(
                                            rownames(df4) %in% mirtarbase[[1]]))
    rownames(df5) <- rownames(corrTable)
    # Which occur in all 3?
    df6 <- df5 %>% mutate(Pred_Fun = Predicted_Interactions + miRTarBase)
    rownames(df6) <- rownames(corrTable)

    # Add to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list("Int_Matrix" = df6)))
    MAE <- c(MAE, MAE2)

return(MAE)
}
