#' @title significantVals
#' @description Filters out genes in each nested dataframe which are not deemed
#' significantly differentially expressed. The results will be stored in
#' the metadata of an MAE object.
#' @param MAE MultiAssayExperiment object. The output of significantVals will
#' be added to this MAE object. It is recommended to use the MAE object
#' created from genesList here.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param geneList A list of nested dataframes if 'c' or a list of lists
#' with nested dataframes if 's'. This will be the output from the
#' genesList function. The resulting list will be found as metadata of
#' the MAE used in the genesList function.
#' @param maxVal Numeric value which represents the maximum cut off value
#' for significance e.g. 0.05.
#' @param stringVal Character. Common DE result type in all nested
#' dataframes which should be used for filtration e.g. pval, adjPval,
#' qval. Make sure this matches the colnames.
#' @return A list with only significantly differentially expressed genes
#' which will also be stores as metadata of a MAE.
#' @export
#' @usage significantVals(MAE, method = '', geneList, maxVal, stringVal = '')
#' @examples
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' Data <- startObject(miR = miR, mRNA = mRNA)
#'
#' Data <- combineGenes(MAE = Data, miR_data = assay(Data, 1),
#'                      mRNA_data = assay(Data, 2))
#'
#' Data <- genesList(MAE = Data, method = 'c', genetic_data = assay(Data, 3),
#'                   timeString = 'D')
#'
#' Data <- significantVals(MAE = Data, method = 'c',
#'                         geneList = metadata(Data)[[1]],
#'                         maxVal = 0.05, stringVal = "adjPVal")
significantVals <- function(MAE, method, geneList, maxVal, stringVal){

    if (missing(method)) stop('Add a MultiAssayExperiment. Results will
                              be added to the MAE object. Please use the
                              genesList function before significantVals.')

    if (missing(method)) stop('method should be "s" for separate analysis and
                              "c" for combined analysis.')

    if(missing(geneList)) stop('Input list of miR and mRNA data. Please use
                               the genesList function first. Results of
                               genesList should be stored in the metadata of
                               the MAE used in the genesList function.')

    if(missing(maxVal)) stop('Add number as cutoff threshold e.g. 0.05.')

    if(missing(stringVal)) stop('Add differential expression result type
                                 to use as filtration point e.g. "log2FC",
                                 "adjPval", "qVal".')

    metadata <- `metadata<-` <- NULL

    if (method == 'c') {

        # Filter for maxVal within the stringVal in dataframes in a list
        filtered_genelist <- lapply(geneList, function(df) df[df[[grep(
            stringVal, names(df),value = TRUE)]] < maxVal,])

        metadata(MAE)[["filtered_genelist"]] <- filtered_genelist

        return(MAE)

    } else if (method == 's') {

        # Filter for maxVal within the stringVal in dataframes in a list of
        # lists
        filtered_genelist <- lapply(geneList, function(ls){
            lapply(ls, function(df) df[df[[grep(stringVal, names(df),
                                                value = TRUE)]] < maxVal,])})

        metadata(MAE)[["filtered_genelist"]] <- filtered_genelist

        return(MAE)

    } else {stop('Please insert method c for combined analysis or s for
                 separate analysis')}
}
