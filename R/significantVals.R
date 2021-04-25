#' @title significantVals
#' @description Filters out genes in each nested dataframe which are not deemed
#' significantly differentially expressed. Each sample will be filtered
#' independently.
#' @param MAE MultiAssayExperiment to store the output of significantVals. It is
#' recommended to use the MAE used in the genesList.
#' @param method Either "c" or "s", respectively for combined or separated
#' analysis.
#' @param geneList A list of nested dataframes if "c" analysis is used or a list
#' of lists of nested dataframes if "s" is used. This will be the output from
#' them genesList function. The resulting list will be found as metadata, in
#' the MAE used in the genesList function.
#' @param maxVal Numeric value which represents the maximum cut off value
#' for significance e.g. 0.05.
#' @param stringVal Character. Common DE result type which is found in all
#' nested dataframes. This will be used for filtration e.g. pval, adjPval or
#' qval. Make sure the spelling matches the colnames for each sample.
#' @return A list of dataframes with only significantly differentially
#' expressed genes. Output will be stored as metadata within the input MAE.
#' @export
#' @usage significantVals(MAE, method = '', geneList, maxVal, stringVal = '')
#' @examples
#'
#' miR <- mm_miR[1:50,]
#'
#' mRNA <- mm_mRNA[1:100,]
#'
#' Data <- startObject(miR = mm_miR, mRNA = mm_mRNA)
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

    if (missing(MAE)) stop('MAE is missing. Add a MultiAssayExperiment. Output of significantVals will be stored in the MAE. Please use the genesList function first')

    if (missing(method)) stop('method is missing. Please add method "c" for combined analysis or "s" for separated analysis')

    if(missing(geneList)) stop('geneList is missing. Add list of nested dataframes. Please use the genesList function first. Results of genesList should be stored in the metadata of the MAE used in the genesList function.')

    if(missing(maxVal)) stop('maxVal is missing. Add number as cutoff threshold e.g. 0.05.')

    if(missing(stringVal)) stop('stringVal is missing. Add differential expression result type to use for filtration e.g. "pval", "adjPval", "qVal".')

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

    } else print('Enter c or s as method.')
}
