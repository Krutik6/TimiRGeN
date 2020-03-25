#' @title AddPrefix
#' @description Adds 'miR_' or 'mRNA_' prefix to colnames of dataframes. Can
#' also add any other prefix, if user has other datatypes to explore. Colnames
#'should end up in the following naming type: 'genetype_timepoint.datatype'.
#' @param MAE MultiAssayExperiment
#' @param gene_df mRNA or miR results from differential expression analysis 
#' inside of a MAE.
#' @param prefixString miR_ or mRNA_
#' @return mRNA or miR results now should have mRNA_ or miR_ at the start
#' of each column name
#' @export
#' @usage AddPrefix(MAE, gene_df, prefixString = '')
#' @examples
#' library(MultiAssayExperiment)
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' 
#' MAE <- StartObject(miR = miR, mRNA = mRNA)
#' 
#' MAE <- AddPrefix(MAE = MAE, gene_df = assay(MAE, 1), 
#'                  prefixString = "miR")
#' 
#' MAE <- AddPrefix(MAE = MAE, gene_df = assay(MAE, 2), 
#'                  prefixString = "mRNA")
AddPrefix <- function(MAE, gene_df, prefixString){
    if (missing(MAE)) stop('Add MultiAssayExperiment');
    if (missing(gene_df)) stop('Input as.data.frame: rows as genenames and 
                                columns as timepoints and differential 
                                expression results. Remember colnames should be 
                                in this format: timepoint.resulttype e.g.
                                D1.log2FC');
    if (missing(prefixString)) stop('Input a genetpye prefix e.g. miR or mRNA');

    gene_df <- as.data.frame(gene_df) 
    # If gene_df already has the prefixString as a prefix then == NO
    ifelse(test = grepl(prefixString, names(gene_df)) == FALSE,
            yes = colnames(gene_df) <- paste(prefixString, 
                                            colnames(gene_df), sep = '_'),
            no = print('miR/mRNA info is fine'))
    # Add to MAE
    MAE2 <- suppressMessages(MultiAssayExperiment(experiments = list(
                                                                x = gene_df)))
    
    x <- paste(prefixString, "p", sep = "_")
    names(MAE2) <- x
    
    MAE <- c(MAE, MAE2)

return(MAE)
}
