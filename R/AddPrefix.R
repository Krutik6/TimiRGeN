#' @title AddPrefix
#' @description Adds 'miR_' or 'mRNA_' prefix to colnames of dataframes. Can
#' also add any other prefix, if user has other datatypes to explore. Colnames
#'should end up in the following naming type: 'genetype_timepoint.datatype'.
#' @param gene_df mRNA or miR results from differential expression analysis 
#' inside of a MAE.
#' @param prefixString miR_ or mRNA_
#' @return mRNA or miR results now should have mRNA_ or miR_ at the start
#' of each column name
#' @export
#' @usage AddPrefix(gene_df, prefixString = '')
#' @examples
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' AddPrefix(gene_df = MAE@ExperimentList$miR,
#' prefixString = "miR" ) -> MAE@ExperimentList$miR_p
#' AddPrefix(gene_df = MAE@ExperimentList$mRNA, 
#' prefixString = "mRNA") -> MAE@ExperimentList$mRNA_p
AddPrefix <- function(gene_df, prefixString){
if (missing(gene_df)) stop('Input as.data.frame: rows as genenames and
columns as timepoints and differential expression results. Remember colnames 
should be in this format: timepoint.resulttype e.g. D1.log2FC');
if (missing(prefixString)) stop('Input a genetpye prefix e.g. miR or mRNA');
gene_df <- as.data.frame(gene_df) 
ifelse(test = grepl(prefixString, names(gene_df)) == FALSE,
yes = colnames(gene_df) <- paste(prefixString, colnames(gene_df), sep = '_'),
no = print('miR/mRNA info is fine'))
return(as.data.frame(gene_df))
}
