#' @title DataMiningMatrix
#' @description Mines out predicted/ functional interactions which correspond
#'between Targetscans, miRDB, miRTarBase and the interactions data frame.
#' @param interactionsDF Correlation matrix between the mRNAs from the
#' wikipathway of interest and miRNA data.
#' @param targetscanInt Targetscans predicted interactions.
#' @param mirdbInt miRDB predcited interactions.
#' @param mirtarbaseInt miRTarBase functional interactions.
#' @return A matrix of information which cross references the occurance of
#'microRNA-mRNA interactions between databases and the given data.
#' @export
#' @importFrom stringr %>%
#' @import dplyr
#' @importFrom dplyr mutate
#' @usage DataMiningMatrix(interactionsDF, targetscanInt ,
#'mirdbInt, mirtarbaseInt)
DataMiningMatrix <- function(interactionsDF, targetscanInt, mirdbInt,
        mirtarbaseInt){
        if (missing(interactionsDF)) stop('Input interactions
        (miR-mRNA(wikipathway)) dataframe.');
        if (missing(targetscanInt)) stop('Input filtered targetscans file,
        interactions column.');
        if (missing(mirdbInt)) stop('Input filtered targetscans file,
        interactions column.');
        if (missing(mirtarbaseInt)) stop('Input filtered targetscans file,
        interactions column.');
        TargetScan <- miRDB <- Predicted_Interactions <- miRTarBase <- NULL
        interactionsDF %>% mutate(TargetScan = as.integer(
        rownames(interactionsDF) %in% targetscanInt)) -> df2
        rownames(df2) <- rownames(interactionsDF)
        df2 %>% mutate(miRDB = as.integer(rownames(df2) %in% mirdbInt)) -> df3
        rownames(df3) <- rownames(interactionsDF)
        df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB) -> df4
        rownames(df4) <- rownames(interactionsDF)
        df4 %>% mutate(miRTarBase = as.integer(
        rownames(df4) %in% mirtarbaseInt)) -> df5
        rownames(df5) <- rownames(interactionsDF)
        df5 %>% mutate(Pred_Fun = Predicted_Interactions + miRTarBase) -> df6
        rownames(df6) <- rownames(interactionsDF)
return(df6)
}
