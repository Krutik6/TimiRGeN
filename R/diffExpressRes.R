#' @title diffExpressRes
#' @description diffExpressRes will produce a dataframe which contains data
#' for only one result type, along with an ID of choice. It is recommended to
#' use this function on a DE results which represents abundance such as log2fc
#' or average expression, as this data will be averaged and correlated later in
#' the analysis. This is to be used for miR and mRNA data individually.
#' @param MAE MultiAssayExperiment to store the output of diffExpressRes
#' within it. This function is to be used after pathways of interest have been
#' identified by enrichWiki or returnCluster. It is recommended to
#' store all diffExpressRes results in the MAE used in enrichWiki and/ or
#' returnCluster.
#' @param df mRNA or miR dataframe (rownames as genes and DE results as
#' columns). These will be found as assays in the MAE object used within the
#' startObject function.
#' @param dataType Column name to take an average from e.g. "Log2FC", "AveExp".
#'  This string should be found consistently in the column names of your input
#'  data. It is recommended to use a DE result value which represents abundance,
#'  rather than confidence.
#' @param genes_ID Dataframe that was created from a getIds function e.g.
#' mRNA_ensembl or miR_entrez. Use the same ID type for miR and mRNA data.
#' These dataframes will be found as assays within the MAE which stores results
#' from the getIds functions.
#' @param idColumn Name of column to use as the merge point. If Column names in
#' getIds results have not been changed, it should be "GENENAME". Default has
#' been left as "GENENAME".
#' @param name New name of the assay. Should be a unique string. Remember
#' each assay in a MAE must have a unique name.
#' @return Dataframe with only a single result type from DE (e.g. Log2FC) and
#' an ID type e.g. entrezIDs. Output will be stored as an assay in the input
#' MAE.
#' @export
#' @usage diffExpressRes(MAE, df, dataType = '', genes_ID, idColumn = '',
#' name = '')
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR[1:100,]
#'
#' mRNA <- mm_mRNA[1:200,]
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' MAE <- getIdsMrna(MAE, assay(MAE, 2), "useast", 'mmusculus')
#'
#' MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
#'                      genes_ID = assay(MAE, 7),
#'                      idColumn = 'GENENAME',
#'                      name = "mRNA_log2fc")
diffExpressRes <- function(MAE, df, dataType, genes_ID, idColumn = 'GENENAME',
                           name){

    if (missing(MAE)) stop('MAE is missing. Add a MAE. This will store output from diffExpressRes. Please use the startObject and getIds functions first.')

    if (missing(df)) stop('df is missing. Add miR or mRNA dataframe which contains genes and results from DE. Please use the startObject function before diffExpressRes. Results of startObject will be found as assays in the MAE used in the startObject function.')

    if (missing(dataType)) stop('dataType is missing. Add string to represent a common result from DE e.g "AveExp", "Log2fc".')

    if (missing(genes_ID)) stop('genes_ID is missing. Add dataframe from getIDs functions. One column name is named "GENENAME", the other contain entrez or ensembl IDs. Please use getIdsMir and getIdsMrna functions first. Output of getIds functions will be stored as assays in the MAE used in the getIds function.')

    if (missing(name)) stop('name is missing. Add name of new data frame. This should be a unique string.')

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

    # Change names of new MAE2 object
    names(MAE2) <- name

    MAE <- c(MAE, MAE2)

    return(MAE)
}
