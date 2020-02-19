#' @title SignificantVals
#' @description Filter out genes in each nested dataframe which are not deemed
#' significantly differentially expressed.
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
#' @usage SignificantVals(method = '', geneList, maxVal, stringVal = '')
#' @examples
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> Data
#' 
#' CombineGenes(miR_data = Data@ExperimentList$miR, mRNA_data = 
#' Data@ExperimentList$mRNA) -> Data@ExperimentList$genetic_data
#' 
#' GenesList(method = 'c', genetic_data = Data@ExperimentList$genetic_data,
#' timeString = 'D') -> Data@metadata$genelist
#' 
#' SignificantVals(method = "c", geneList = Data@metadata$genelist, 
#' maxVal = 0.05, stringVal = "adjPVal") -> Data@metadata$filtered_genelist
SignificantVals <- function(method, geneList, maxVal, stringVal){
if (missing(method)) stop('method should be s for separate analysis and
c for combined analysis.')
if (method == 'c') {
if(missing(geneList)) stop('Input list of miR and mRNA data.');
if(missing(maxVal)) stop('Input integer as cutoff threshold e.g. 0.05.');
if(missing(stringVal)) stop('Input differential expression result type
to use as filtration point e.g. log2FC,
adjPval, qVal.');
X <- lapply(geneList, function(df) df[df[[grep(stringVal, names(df),
value = TRUE)]] < maxVal,])
return(X)
} else if (method == 's') {
if(missing(geneList)) stop('Input list of miR and mRNA data.');
if(missing(maxVal)) stop('Input integer as cutoff threshold e.g. 0.05.');
if(missing(stringVal)) stop('Input differential expression result type
to use as filtration point e.g. log2FC,
adjPval, qVal.');
lapply(geneList, function(ls){
lapply(ls, function(df) df[df[[grep(stringVal, names(df),
value = TRUE)]] < maxVal, ])
})
} else {stop('Please insert method c for combined
analysis or s for seperate analysis')}
}
