#' @title eNames
#' @description Extract the gene IDs from the nested dataframes within the
#' metadata of the MAE object, which is the output of the addIds function.
#' Results of eNames will be stored in the metadata of the resulting
#' MAE object.
#' @param MAE MultiAssayExperiment to store output of eNames.
#' It is recommended to use the MAE object which contains addIds results.
#' @param method Either 'c' or 's' for combined or separated analysis.
#' @param gene_IDs List of DE data and associated gene IDs.
#' Output from addIds function, this should be found in the metadata of the MAE
#' used in the addIds function.
#' @param ID_Column The entrez/ ensembl ID column in each dataframe.
#' Should be the last column. This should be 2+ the number of DE results per
#' time point. e.g. if a user has log2f and adj.P.val results for each
#' time point, then the fourth column will contain the gene ID information.
#' @return A single list of entrez/ ensembl IDs for each dataframe
#'within the lists.
#' @export
#' @importFrom stats complete.cases
#' @usage eNames(MAE, method = '', gene_IDs, ID_Column)
#' @examples
#'library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' Data <- startObject(miR = miR, mRNA = mRNA)
#'
#' Data <- getIdsMirMouse(Data, assay(Data, 1))
#'
#' Data <- getIdsMrnaMouse(Data, assay(Data, 2), mirror = 'www')
#'
#' Data <- combineGenes(MAE = Data, miR_data = assay(Data, 1),
#'                      mRNA_data = assay(Data, 2))
#'
#' Data <- genesList(MAE = Data, method = 'c', genetic_data = assay(Data, 9),
#'                   timeString = 'D')
#'
#' Data <- significantVals(MAE = Data, method = 'c',
#'                         geneList = metadata(Data)[[1]],
#'                         maxVal = 0.05, stringVal = "adjPVal")
#'
#' Data <- addIds(MAE = Data, method = "c",
#'               filtered_genelist = metadata(Data)[[2]],
#'               miR_IDs = assay(Data, 3), mRNA_IDs = assay(Data, 7))
#'
#' Data <- eNames(MAE = Data, method = "c", gene_IDs = metadata(Data)[[3]],
#'                ID_Column = 4)
eNames <- function(MAE, method, gene_IDs, ID_Column){

    if (missing(MAE)) stop('Add MultiAssayExperiment. The results of eNames
                           will be stored in the MAE. Please use the addIds
                           function first.')

    if (missing(method)) stop('method should be "s" for separate analysis and
                               "c" for combined analysis.')

    if (missing(gene_IDs)) stop('Input a list of nested dataframes with
                                 entrezID/ ensembl gene name information.
                                 Please use the addIDs function first. Output
                                 of the addIds function will be stored as
                                 metadata of the MAE used in the addIds
                                function.')

    if (missing(ID_Column)) stop('Input an integer representing the
                                  column which contains gene ID information.
                                  This should be 2+ the number of DE results
                                  per time point.')

    metadata <- `metadata<-` <- NULL

    # If c is choses
    if (method == 'c') {
        # Retreive all the string from the ID_Column from each list
        e <- lapply(gene_IDs, function(x){x[[ID_Column]]})

        ID_list <- lapply(e, function(x){ x[stats::complete.cases(x)]})

        # store into MAE object
        metadata(MAE)[["ID_list"]] <- ID_list

    return(MAE)

    } else if (method == 's') {

        # Retrieve all the string from the ID_Column from each list within
        # a list
        Y <- sapply(gene_IDs, function(x){sapply(x, `[[`, ID_Column)})

        X <- vapply(gene_IDs, function(x){list(names(x))}, FUN.VALUE = list(1))

        Xnames <- unlist(X)

        names(Y) <- Xnames

        ID_list <- lapply(Y, function(x){ x[stats::complete.cases(x)]})

        # store into MAE object
        metadata(MAE)[["ID_list"]] <- ID_list

    return(MAE)

} else {stop('Please insert method c for combined analysis or s for
              separate analysis')}
}
