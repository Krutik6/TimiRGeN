#' @title addPrefix
#' @description Adds 'miR_' or 'mRNA_' prefix to colnames of dataframes. Can
#' also add any other prefix, if user has other datatypes to explore. Colnames
#' should end up in the following naming type: 'genetype_timepoint.datatype'.
#' @param MAE MultiAssayExperiment
#' @param gene_df mRNA or miR results from differential expression analysis
#' inside of a MAE.
#' @param prefixString miR_ or mRNA_
#' @return mRNA or miR results now should have mRNA_ or miR_ at the start
#' of each column name.
#' @export
#' @usage addPrefix(MAE, gene_df, prefixString = '')
#' @examples
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- addPrefix(MAE = MAE, gene_df = assay(MAE, 1),
#'                  prefixString = "miR")
#'
#' MAE <- addPrefix(MAE = MAE, gene_df = assay(MAE, 2),
#'                  prefixString = "mRNA")
addPrefix <- function(MAE, gene_df, prefixString){
    if (missing(MAE)) stop('Add MultiAssayExperiment.')

    if (missing(gene_df)) stop('Input output from startObject.')

    if (missing(prefixString)) stop('Input a genetpye prefix e.g. miR or mRNA.')

    gene_df <- as.data.frame(gene_df)

    # If gene_df already has the prefixString as a prefix then == NO
    ifelse(test = grepl(prefixString, names(gene_df)) == FALSE,
            yes = colnames(gene_df) <- paste(prefixString,
                                             colnames(gene_df),
                                             sep = '_'),
            no = print('miR/mRNA info is fine'))
    # Add to MAE
    MAE2 <- suppressMessages(MultiAssayExperiment(experiments = list(
                                                                x = gene_df)))

    # Add prefix to column names
    x <- paste(prefixString, "p", sep = "_")

    names(MAE2) <- x

    MAE <- c(MAE, MAE2)

return(MAE)
}
