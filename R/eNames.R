#' @title eNames
#' @description Extract the ID names from the nested dataframes within the
#' metadata of the MAE object.
#' @param method Respectively either 'c' or 's' for combined or separated
#'analysis.
#' @param gene_IDs A nested list of dataframes(c)/ list of lists of nested
#'dataframes (s), and each dataframe has a column for entrez or ensembl ID.
#' @param ID_Column The entrez/ ensembl ID column in each dataframe.
#'Should be the last column.
#' @return A single list of entrez/ ensembl IDs for each dataframe
#'within the lists.
#' @export
#' @usage eNames(method = '', gene_IDs, ID_Column)
#' @examples
#'library(biomaRt)
#'library(clusterProfiler)
#'library(org.Mm.eg.db)
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'getIDs_miR_mouse(MAE = MAE, miR = MAE@ExperimentList$miR) -> MAE
#'getIDs_mRNA_mouse(MAE = MAE, mRNA = MAE@ExperimentList$mRNA) -> MAE
#'
#'CombineGenes(miR_data = MAE@ExperimentList$miR,
#'mRNA_data = MAE@ExperimentList$mRNA) -> MAE@ExperimentList$genetic_data
#' 
#'GenesList(method = 'c', genetic_data = MAE@ExperimentList$genetic_data,
#'timeString = 'D') -> MAE@metadata$genelist
#' 
#'SignificantVals(method = "c", geneList = MAE@metadata$genelist,
#'maxVal = 0.05, stringVal = "adjPVal") -> MAE@metadata$filtered_genelist
#' 
#'AddIDs(method = "c", filtered_genelist = MAE@metadata$filtered_genelist,
#'miR_IDs = MAE@ExperimentList$miR_entrez,
#'mRNA_IDs = MAE@ExperimentList$mRNA_entrez) -> MAE@metadata$gene_entrez
#' 
#'eNames(method = "c", gene_IDs = MAE@metadata$gene_entrez, ID_Column = 4
#') -> MAE@metadata$e_list
eNames <- function(method, gene_IDs, ID_Column){
if (missing(method)) stop('method should be s for separate analysis and
c for combined analysis.')
if (missing(gene_IDs)) stop('Input a list of nested dataframes with
entrezID/ ensembl gene name information.');
if (missing(ID_Column)) stop('Input an interger representing the
column which contains entrezIDs
information');
if (method == 'c') {
lapply(gene_IDs, function(x){x[[ID_Column]]}) -> e
lapply(e, function(x){ x[complete.cases(x)]}) -> y
return(y)
} else if (method == 's') {
sapply(gene_IDs, function(x){sapply(x, `[[`, ID_Column)}) -> Y
vapply(gene_IDs, function(x){list(names(x))}, FUN.VALUE = list(1)) -> X
unlist(X) -> Xnames
names(Y) <- Xnames
lapply(Y, function(x){ x[complete.cases(x)]}) -> y
return(y)
} else {stop('Please insert method c for combined analysis or s for
seperate analysis')}
}
