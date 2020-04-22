#' @title createClusters
#' @aliases createClusters
#' @description Create soft clusters to assess how your datasets change in
#' gene abundance during your time course in different wikipathways.
#' This function is to be used before ClusterCheck and Quickfuzz.
#' @param MAE MultiAssayObject where the cluster information will be stored.
#' @param method Either C or S for combined or separate analysis.
#' @param percentMatrix Output after TurnPercent.
#' @param noClusters Number of clusters to create, the default is 5.
#' @param dataString Only for use in S analysis. Insert the prefix string e.g.
#' mRNA or miR.
#' @param variance Numeric from 0-1 to control strictness of filtering. Higher
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
#' createClusters(MAE, method, percentMatrix, noClusters, dataString = '',
#'                variance)
#' @examples
#' library(MultiAssayExperiment)
#' library(Mfuzz)
#' MAE <- MultiAssayExperiment()
#' metadata(MAE)[["e_list"]] <- e_list
#' metadata(MAE)[["w_list"]] <- w_list[1:10]
#' MAE <- wikiMatrix(MAE, ID_list = metadata(MAE)[[1]],
#'                   wp_list = metadata(MAE)[[2]])
#'
#' MAE <- turnPercent(MAE = MAE,
#'                    wikiMatrix = assay(MAE, 1),
#'                    rowInt = 6)
#'
#'MAE <- createClusters(MAE, method = "c",
#'                     percentMatrix = assay(MAE, 2),
#'                     noClusters = 2, variance = 0.99)
createClusters <- function(MAE, method, percentMatrix, noClusters = 5,
                           dataString, variance = 0){

    if (missing(MAE)) stop('Add MultiAssayExperiment');
    if (missing(method)) stop('Enter c for combined or s for separate');
    if (missing(percentMatrix)) stop('Dataframe which is the output from
                                      TurnPercent function');

    metadata <- `metadata<-` <- NULL

    # ready the percentMatrix
    X <- percentMatrix
    df <- as.data.frame(t(X))
    df$Total <- NULL
    # Convert factors into numeric
    df <- data.matrix(frame = df, rownames.force = NA)

    df <- round(df, 0)
    df <- na.omit(df)
    # If == s then subset data by common string e.g mRNA, miR
    if (method == 's') {
    df2 <- df[, grepl(dataString, names(df))]
    # If == c do nothing
    } else if (method == 'c') {
    df2 <- df

    } else print('Select s for separated analysis or c for combined
    analysis')
    # standardise data

    Eset <- new('ExpressionSet', exprs = as.matrix(df2))

    Eset_sd <- filter.std(Eset, min.std = variance)

    Eset_st <- standardise(Eset_sd)

    m <- mestimate(Eset_st)

    # Perfrom mfuzz
    Clusters <- mfuzz(Eset_st, centers = noClusters, m=m)

    X <- as.data.frame(Clusters$membership)

    # retreive pathway names
    for (i in seq_along(rownames(X))) {
        getPathwayInfo(rownames(X)[i])[[3]]
        } ->  X$Description[i]

    # Save data in MAE
    MAE2 <- suppressMessages(MultiAssayExperiment(list("ClusterData" = X,
                                                       "MfuzzData" = Eset_st)))
    metadata(MAE2)[["Clusters"]] <- Clusters
    MAE <- c(MAE, MAE2)

return(MAE)
}
