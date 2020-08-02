#' @title createClusters
#' @description Create soft clusters to assess change in gene abundance during
#' the time course in different wikipathways. createClusters will create
#' 3 data files. Clusters which will be stored in the metadata,
#' MfuzzData which contains fuzzy clustering information will be stored
#' as an experiment and ClusterData which gives cluster-pathway fit information
#' which is stored as an assay. This function may take some time as it
#' downloads pathway information.
#' @param MAE MultiAssayExperiment where the results from createClusters will be
#' stored. It is advised to use the MAE object which stores the output of
#' by turnPercent.
#' @param method Either "c" or "s" for combined or separate analysis.
#' @param percentMatrix A matrix containing wikipathway-data information. It
#' is output from the turnPercent function and will be stored as an assay
#' within the MAE used in the turnPercent function.
#' @param noClusters Number of clusters to create, the default is 5.
#' @param dataString Only for use in "s" analysis. Insert the prefix string e.g.
#' mRNA or miR. The string added should be the same prefix added using the
#' addPrefix function.
#' @param variance Numeric value from 0-1 to control strictness of filtering.
#' Higher variance means more pathways will be excluded from the analysis.
#' @return Clusters(metadata): A list to be used as the input in checkClusters.
#' MfuzzData(ExperimentList): An ExpressionSet object to be input for
#' quickFuzz.
#' ClusterData(assay): A breakdown of how each pathway fitted with
#' each cluster and is the input for returnCluster.
#' @export
#' @importFrom Mfuzz filter.std standardise mestimate mfuzz
#' @importFrom stats na.omit
#' @importFrom rWikiPathways getPathwayInfo
#' @usage createClusters(MAE, method, percentMatrix, noClusters,
#' dataString = '', variance)
#' @examples
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list_mouse
#'
#' metadata(MAE)[["w_list"]] <- w_list_mouse[1:10]
#'
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

    if (missing(MAE)) stop('Add MultiAssayExperiment. Data from createClusters
                           will be stored in this MAE object. Please use
                           turnPercent first.')

    if (missing(method)) stop('Enter "c" for combined or "s" for separate')

    if (missing(percentMatrix)) stop('Dataframe which contains pathways-samples
                                     information as percentages. Please use the
                                     turnPercent function first. Results of this
                                     will be stored as an assay within the MAE
                                     used in the turnPercent function.')

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

      if (missing(dataString)) stop('Add prefix e.g. "mRNA" or "miR".')

      df <- as.data.frame(df)

      df2 <- df[, grepl(dataString, names(df))]

    # If == c do nothing
    } else if (method == 'c') {

      df2 <- df

    } else print('Select s for separated analysis or c for combined
    analysis')

    # standardise data using standard Mfuzz code
    Eset <- new('ExpressionSet', exprs = as.matrix(df2))

    Eset_sd <- Mfuzz::filter.std(Eset, min.std = variance)

    Eset_st <- Mfuzz::standardise(Eset_sd)

    m <- Mfuzz::mestimate(Eset_st)

    # Perform mfuzz
    Clusters <- Mfuzz::mfuzz(Eset_st, centers = noClusters, m=m)

    # Get cluster membership information
    X <- as.data.frame(Clusters$membership)

    # retrieve pathway names using rWikipathways using membership
    for (i in seq_along(rownames(X))) {
      X$Description[i] <- rWikiPathways::getPathwayInfo(rownames(X)[i])[[3]]
        }

    # Save data in MAE
    MAE2 <- suppressMessages(MultiAssayExperiment(list("ClusterData" = X,
                                                       "MfuzzData" = Eset_st)))
    metadata(MAE2)[["Clusters"]] <- Clusters

    MAE <- c(MAE, MAE2)

return(MAE)
}
