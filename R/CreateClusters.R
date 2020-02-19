#' @title CreateClusters
#' @aliases CreateClusters
#' @description Create soft clusters to assess how your datasets change in
#' gene abundance during your time course in different wikipathways.
#' This function is to be used before ClusterCheck and Quickfuzz.
#' @param method Either C or S for combined or separate analysis.
#' @param MAE MultiAssayObject where the cluster information will be stored.
#' @param Percent_matrix Output after TurnPercent.
#' @param no.clusters Number of clusters to create, the default is 5.
#' @param Data_string Only for use in S analysis. Insert the prefix string e.g.
#' mRNA or miR.
#' @param Variance Numeric from 0-1 to control strictness of filtering. Higher
#' Variance means wikipathways that are undergoing a higher change in the
#' data will be clustered.
#' @importFrom Mfuzz filter.std standardise mestimate mfuzz
#' @return Clusters(metadata): A list to be used as the input in plot functions.
#' Mfuzzdata: An input for QuickFuzz. Mfuzzdata(ExperimentList): An
#' ExpressionSet objectto be input for Mfuzz clustering. 
#' ClusterData(ExperimentList): A breakdown of how each pathway fitted with 
#' each cluster.
#' @export
#' @usage 
#' CreateClusters(method, MAE, Percent_matrix, no.clusters, Data_string = '',
#' Variance)
#' @examples
#'library(Mfuzz)
#'library(rWikiPathways)
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'e_list -> MAE@metadata$elist 
#'w_list[1:10] -> MAE@metadata$wlist
#'WikiMatrix(e_list = MAE@metadata$elist , 
#'wp_list = MAE@metadata$wlist) -> MAE@ExperimentList$Wmat
#'
#'TurnPercent(wikiMatrix = MAE@ExperimentList$Wmat,
#'rowInt = 4) -> MAE@ExperimentList$Pmat
#'
#'CreateClusters(method = "c", MAE, Percent_matrix = MAE@ExperimentList$Pmat,
#' no.clusters = 2, Variance = 0.99) -> MAE
CreateClusters <- function(method, MAE, Percent_matrix, no.clusters = 5,
Data_string, Variance = 0){
as.data.frame(t(Percent_matrix)) -> df
df$Total <- NULL
df[vapply(df, is.factor, logical(1))] <- lapply(df[vapply(
df, is.factor, logical(1))], function(x) as.numeric(as.character(x)))
round(df, 0) -> df
na.omit(df) -> df
if (method == 's') {
df[, grepl(Data_string, names(df))] -> df2
} else if (method == 'c') {
df -> df2
} else print('Select s for separated analysis or c for combined
analysis')
Eset <- new('ExpressionSet', exprs = as.matrix(df2))
Eset_sd <- filter.std(Eset, min.std = Variance)
Eset_st <- standardise(Eset_sd)
m <- mestimate(Eset_st)
cl <- mfuzz(Eset_st, centers = no.clusters, m=m)
MAE@ExperimentList$MfuzzData <- Eset_st
MAE@metadata$Clusters <- cl
as.data.frame(cl$membership) -> X
for (i in seq_along(rownames(X))) {
getPathwayInfo(rownames(X)[i])[[3]]
} -> X$Description[i]
MAE@ExperimentList$ClusterData <- X
return(MAE)
}
