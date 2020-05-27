#' @title eNames
#' @description Extract the ID names from the nested dataframes within the
#' metadata of the MAE object.
#' @param MAE MultiAssayExperiment Object.
#' @param method Either 'c' or 's' for combined or separated analysis.
#' @param gene_IDs Output from addIDs function.
#' @param ID_Column The entrez/ ensembl ID column in each dataframe.
#'Should be the last column.
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
#' Data <- combineGenes(MAE = Data, miR_data = assay(Data, 1),
#'                      mRNA_data = assay(Data, 2))
#'
#' Data <- genesList(MAE = Data, method = 'c', genetic_data = assay(Data, 3),
#'                   timeString = 'D')
#'
#' Data <- significantVals(MAE = Data, method = 'c',
#'                         geneList = metadata(Data)[[1]],
#'                         maxVal = 0.05, stringVal = "adjPVal")
#'
#' Data <- addIds(MAE = Data, method = "c",
#'               filtered_genelist = metadata(Data)[[2]],
#'               miR_IDs = assay(Data, 3), mRNA_IDs = assay(Data, 3))
#'
#' Data <- eNames(MAE = Data, method = "c", gene_IDs = metadata(Data)[[3]],
#'                ID_Column = 4)
eNames <- function(MAE, method, gene_IDs, ID_Column){

    if (missing(MAE)) stop('Add MultiAssayExperiment.')

    if (missing(method)) stop('method should be s for separate analysis and
                               c for combined analysis.')

    if (missing(gene_IDs)) stop('Input a list of nested dataframes with
                                 entrezID/ ensembl gene name information.')

    if (missing(ID_Column)) stop('Input an interger representing the
                                  column which contains entrezIDs information.')

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

        # Retreive all the string from the ID_Column from each list within
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
seperate analysis')}
}
