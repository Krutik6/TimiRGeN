#' @title significantVals
#' @description Filter out genes in each nested dataframe which are not deemed
#' significantly differentially expressed.
#' @param MAE MultiAssayExperiment object.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param geneList A list of nested dataframes if 'c' or a list of lists
#' with nested dataframes if 's'.
#' @param maxVal Integer. Rows will be kept only if they are
#' lower than this value.
#' @param stringVal Character. Common resulttype in all nested
#' dataframes which should be point of filtration e.g. pval, adjPval,
#' qval. Make sure this matches the colnames.
#' @return A list in a similair structure but with only significantly
#' differentially expressed genes which will also be stores in the
#' metadata area of an MAE object.
#' @export
#' @usage significantVals(MAE, method = '', geneList, maxVal, stringVal = '')
#' @examples
#' miR <- mm_miR
#' mRNA <- mm_mRNA
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

    if (missing(method)) stop('Use MultiAssayExperiment.');
    if (missing(method)) stop('method should be s for separate analysis and
                              c for combined analysis.');
    if(missing(geneList)) stop('Input list of miR and mRNA data.');
    if(missing(maxVal)) stop('Input integer as cutoff threshold e.g. 0.05.
                                 ');
    if(missing(stringVal)) stop('Input differential expression result type
                                    to use as filtration point e.g. log2FC,
                                    adjPval, qVal.');

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
                 seperate analysis')}
}
