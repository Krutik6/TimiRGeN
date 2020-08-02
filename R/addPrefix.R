#' @title addPrefix
#' @description Adds 'miR_' or 'mRNA_' prefix to colnames of dataframes. Can
#' also add any other prefix, if user has other data types to explore. Colnames
#' should end up in the following naming system: 'genetype_timepoint.datatype'.
#' This function is essential for separate analysis of miR-mRNA DE data. If
#' using the combined analysis, there is no need to use addPrefix.
#' @param MAE MultiAssayExperiment to store output of addPrefix.
#' It is recommended to use the MAE object which stores output from startObject.
#' @param gene_df mRNA or miR results from differential expression analysis
#' inside of a MAE. Would be stored as an assay within the MAE used in the
#' startObject function.
#' @param prefixString prefix to be added e.g. "miR" or "mRNA".
#' @return Dataframe which has a specific prefix infront of each column name.
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
    if (missing(MAE)) stop('Add MAE to store output of addPrefix. Please
                           use startObject first. ')

    if (missing(gene_df)) stop('Input dataframe containing miR or mRNA DE
                               data. Please use startObject first. Output of
                               startObject should be stored as metadata within
                               the MAE used in startObject.')

    if (missing(prefixString)) stop('Input a prefix string e.g. "miR" or
                                    "mRNA."')

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
