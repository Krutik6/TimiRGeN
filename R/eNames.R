#' @title eNames
#' @description Extract the ID names from the nested dataframes.
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
#' library(biomaRt)
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' miR <- mm_miR
#' miR <- miR[1:100,]
#' mRNA <- mm_mRNA
#' mRNA <- mRNA[1:200,]
#' getIDs_miR_mouse(miR = miR)
#' getIDs_mRNA_mouse(mRNA = mRNA)
#' CombineGenes(miR_data = miR, mRNA_data = mRNA) -> genetic_data
#' GenesList(method = 'c', genetic_data = genetic_data,
#' timeString = 'D') -> genelist
#' SignificantVals(method = 'c', geneList = genelist, maxVal = 0.05,
#' stringVal = 'adjPVal') -> filtered_genelist
#' AddIDs(method = 'c', filtered_genelist = filtered_genelist,
#' miR_IDs = miR_entrez, mRNA_IDs = mRNA_entrez) -> gene_entrez
#' eNames(method = 'c', gene_IDs = gene_entrez, ID_Column = 4) -> e_list
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
