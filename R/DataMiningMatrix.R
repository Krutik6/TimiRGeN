#' @title DataMiningMatrix
#' @description Mines out predicted/ functional interactions which correspond
#'between Targetscans, miRDB, miRTarBase and the interactions data frame.
#' @param corrTable Correlation matrix between the mRNAs from the
#' wikipathway of interest and miRNA data.
#' @param targetscan Targetscans predicted interactions.
#' @param mirdb miRDB predcited interactions.
#' @param mirtarbase miRTarBase functional interactions.
#' @return A matrix of information which cross references the occurance of
#'microRNA-mRNA interactions between databases and the given data.
#' @export
#' @importFrom stringr %>%
#' @import dplyr
#' @importFrom dplyr mutate
#' @usage DataMiningMatrix(corrTable, targetscan , mirdb, mirtarbase)
DataMiningMatrix <- function(corrTable, targetscan, mirdb,
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
corrTable %>% mutate(TargetScan = as.integer(
rownames(corrTable) %in% targetscan[[1]])) -> df2
rownames(df2) <- rownames(corrTable)
df2 %>% mutate(miRDB = as.integer(rownames(df2) %in% mirdb[[1]])) -> df3
rownames(df3) <- rownames(corrTable)
df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB) -> df4
rownames(df4) <- rownames(corrTable)
df4 %>% mutate(miRTarBase = as.integer(
rownames(df4) %in% mirtarbase[[1]])) -> df5
rownames(df5) <- rownames(corrTable)
df5 %>% mutate(Pred_Fun = Predicted_Interactions + miRTarBase) -> df6
rownames(df6) <- rownames(corrTable)
return(df6)
}
