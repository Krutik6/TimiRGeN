#' @title createClusters2
#' @description Creates clusters from a dataframe of mRNAs and miRs. This
#' function should primarily be used when analysing data that has not gone
#' through pair-wise DE. This function will create clusters from longitudinal
#' temporal patterns. createClusters2 will create 3 data files.
#' 1) Clusters will contain cluster logistics information and will be stored as
#' metadata, 2) MfuzzData will contain fuzzy clustering information and will be
#' stored as an experiment, 3) ClusterData will contain cluster-pathway fit
#' information and will be stored as an assay.
#' @param MAE MultiAssayExperiment which will store the results from
#' createClusters2.
#' @param genetic_data A dataframe with miR and mRNA
#' information together. This is the output from the combineGenes function
#' and will be stored as an assay within the MAE used in the combineGenes
#' function.
#' @param noClusters How many clusters should be generated? Default is 5.
#' @return 3 new objects in the input MAE.
#' Clusters(metadata): A list to be used as the input in checkClusters
#' and quickFuzz.
#' MfuzzData(ExperimentList): An ExpressionSet object to be used as input for
#' quickFuzz.
#' ClusterData(assay): An assay to be used as input for returnCluster.
#' @export
#' @usage createClusters2(MAE, genetic_data, noClusters)
#' @examples
#' data(long_data)
#' miRNA <- long_data[c(1:105),]
#' mRNA <- long_data[-c(1:105),]
#'
#' MAE <- startObject(miR = miRNA, mRNA = mRNA)
#'
#' MAE <- combineGenes(MAE, miR_data = assay(MAE, 1),
#'                     mRNA_data = assay(MAE, 2))
#'
#' MAE <- createClusters2(MAE = MAE, genetic_data = assay(MAE, 3))
createClusters2 <- function(MAE, genetic_data, noClusters = 5){

  if (missing(MAE)) stop('MAE is missing. Add MAE. Please use combineGenes first.')

  if (missing(genetic_data)) stop('genetic_data is missing. Input combined miR and mRNA data. Colnames structure should be timepoint.counttype. Please use the combineGenes function first. Output of combineGenes should be stored as an assay within the MAE used in the combineGenes function.')

  colnames(genetic_data) <- as.integer(gsub(colnames(genetic_data),
                                            pattern = "[^0-9.-]",
                                            replacement = ""))

  colnames(genetic_data) <- as.integer(colnames(genetic_data))

  metadata <- `metadata<-` <- NULL

  df <- data.matrix(frame = genetic_data, rownames.force = NA)

  df <- round(df, 0)

  df <- na.omit(df)

  Eset <- new("ExpressionSet", exprs = as.matrix(df))

  Eset_st <- Mfuzz::standardise(Eset)

  m <- Mfuzz::mestimate(Eset_st)

  Clusters <- Mfuzz::mfuzz(Eset_st, centers = noClusters, m = m)

  cs <- as.data.frame(Clusters$membership)

  MAE2 <- suppressMessages(MultiAssayExperiment(list(ClusterData = cs,
                                                     MfuzzData = Eset_st)))

  metadata(MAE2)[["Clusters"]] <- Clusters

  MAE <- c(MAE, MAE2)

  return(MAE)
}
