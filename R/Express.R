#' @title Express
#' @description Get average data for one result type.
#' @param df mRNA or miR dataframe (rownames as genes and experimental
#' variables as columns.)
#' @param dataType Column name to take an average from e.g. Log2FC, AveExp,
#' pval.
#' @param genes_ID Dataframe with that was created earlier e.g. mRNA_ensembl
#' or miR_entrez.
#' @param idColumn Name of column to use as the merge point.
#' @return Returns the average of whichever data the user was interested e.g.
#' Average expression or average log2fc.
#' @export
#' @usage Express(df, dataType = '', genes_ID, idColumn = '')
#' @examples
#' library(biomaRt)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' getIDs_miR_mouse(MAE, MAE@ExperimentList$miR) -> MAE
#' getIDs_mRNA_mouse(MAE, MAE@ExperimentList$mRNA, "useast") -> MAE
#' Express(df = MAE@ExperimentList$mRNA, dataType = 'Log2FC',
#' genes_ID = MAE@ExperimentList$mRNA_entrez,
#' idColumn = 'GENENAME') -> mRNA_express
Express <- function(df, dataType, genes_ID, idColumn){
if (missing(df)) stop('Input miR or mRNA dataframe.');
if (missing(dataType)) stop('Input common name for expression columns from
data. e.g AveExp.');
if (missing(genes_ID)) stop('Input dataframe from bitr function. One column
is GENENAME, the other is entrez or ensembl.');
if (missing(idColumn)) stop('Input which column name to use for merge point
e.g. GENENAME');
df <- as.data.frame(df)
genes_ID <- as.data.frame(genes_ID)
exp <- cbind(names = rownames(df), df[,grep(dataType, colnames(df))])
merge(exp, genes_ID, by.x = 'names', by.y = idColumn, all = TRUE) -> merged
rownames(merged) <- merged[[1]]
merged[[1]] <- NULL
return(merged)
}
