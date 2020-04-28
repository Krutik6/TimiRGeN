#' @title diffExpressRes
#' @description Get average data for one result type.
#' @param MAE MultiAssayExperiment Object.
#' @param df mRNA or miR dataframe (rownames as genes and experimental
#' variables as columns.)
#' @param dataType Column name to take an average from e.g. Log2FC, AveExp,
#' pval.
#' @param genes_ID Dataframe with that was created earlier e.g. mRNA_ensembl
#' or miR_entrez.
#' @param idColumn Name of column to use as the merge point.
#' @param name = New name of the dataframe.
#' @return Returns the average of whichever data the user was interested e.g.
#' Average diffExpressResion or average log2fc.
#' @export
#' @usage diffExpressRes(MAE, df, dataType = '', genes_ID, idColumn = '', name = '')
#' @examples
#' miR <- mm_miR[1:100,]
#' mRNA <- mm_mRNA[1:200,]
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
#' MAE <- getIdsMrnaMouse(MAE, assay(MAE, 2), "useast")
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                genes_ID = assay(MAE, 7),
#'                idColumn = 'GENENAME',
#'                name = "mRNA_log2fc")
diffExpressRes <- function(MAE, df, dataType, genes_ID, idColumn, name){

    if (missing(MAE)) stop('Check if you are using the correct MAE object.')
    if (missing(df)) stop('Input miR or mRNA dataframe.');
    if (missing(dataType)) stop('Input common name for expression columns from
                                 data. e.g AveExp.');
    if (missing(genes_ID)) stop('Input dataframe from bitr function. One column
                                 is GENENAME, the other is entrez or ensembl.');
    if (missing(idColumn)) stop('Input which column name to use for merge point
                                 e.g. GENENAME');
    if (missing(name)) stop('Add name of new data frame.')

    df <- as.data.frame(df)
    genes_ID <- as.data.frame(genes_ID)

    # Isolate columns that have the same dataType
    exp <- cbind(names = rownames(df), df[,grep(dataType, colnames(df))])
    # Merge them with IDs
    merged <- merge(exp, genes_ID, by.x = 'names',
                    by.y = idColumn, all = TRUE)

    rownames(merged) <- merged[[1]]
    merged[[1]] <- NULL

    # Add to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list(x = merged)))
    names(MAE2) <- name
    MAE <- c(MAE, MAE2)
    return(MAE)
}
