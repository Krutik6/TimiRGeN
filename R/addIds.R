#' @title addIds
#' @description Adds entrez or ensembl IDs to the nested dataframes within a
#' list(c) or list of lists (s). The IDs are created by the getIdsMir and
#' getIdsMrna functions, both are needed for addIds.
#' @param MAE MultiAssayExperiment to store the output of addIds.
#' It is recommended to use the MAE which stores results from
#' significantVals.
#' @param method  Either "c" or "s", respectively for combined or separated
#' analysis.
#' @param filtered_genelist A list of nested dataframes if 'c' or a list of
#' lists with nested dataframes if 's'. This will be found as metadata within
#' the MAE object used in the significantVals function.
#' @param miR_IDs miR_ensembl or miR_entrez. Use a getIDsMir function to
#' acquire this. This will be stored as an assay in the MAE used in a getIdsMir
#' function.
#' @param mRNA_IDs mRNA_ensembl or mRNA_entrez. Use a getIDsMrna function to
#' acquire this. This will be stored as an assay in the MAE used in a getIdsMrna
#' function.
#' @return List of dataframes with entrezIDs/ ensembl IDs and gene names
#' as columns which will be stored as metadata in the input MAE.
#' @export
#' @usage addIds(MAE, method, filtered_genelist, miR_IDs, mRNA_IDs)
#' @examples
#'library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' Data <- startObject(miR = miR, mRNA = mRNA)
#'
#' Data <- getIdsMir(Data, assay(Data, 1), orgDB = org.Mm.eg.db, 'mmu')
#'
#' Data <- getIdsMrna(Data, assay(Data, 2), mirror = 'useast', 'mmusculus')
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
addIds <- function(MAE, method, filtered_genelist, miR_IDs, mRNA_IDs){

    if (missing(MAE)) stop('MAE is missing. Add MAE to store output of addIds. Please use significantVals and getIds functions first.')

    if (missing(method)) stop('method is missing. Please add method "c" for combined analysis or "s" for separated analysis')

    if (missing(filtered_genelist)) stop('filtered_genelist is missing. Add filtered list of genes which is listed by genetype and time (s) or just by time (c). Please use significantVals first. Output of significantVals should be stored as metadata within the MAE used in the significantVals function.')

    if (missing(miR_IDs)) stop('miR_IDs is missing. Add dataframe of miR gene IDs. Please use getIdsMirHuman or getIdsMirMouse first. Output of a getIdsMir function should be stored as assays within the MAE used in the getIdsMir function.')

    if (missing(mRNA_IDs)) stop('mRNA_IDs is missing. Add dataframe of mRNA gene IDs. Please use getIdsMrnaHuman or getIdsMrnaMouse first. Output of getIds function should be stored as an assay within the MAE used in the getIdsMrna function.')

    metadata <- `metadata<-` <- NULL

    # For combined analysis
    if (method == 'c') {

        # Extract the entrez/ ensembl data
        geneIDs <- rbind(miR_IDs, mRNA_IDs)

        genes_id <- geneIDs[! duplicated(geneIDs[[1]]),]

        # For each list add the entrez/ ensembl IDs
        X <- lapply(filtered_genelist,
                    function(x){cbind('GENENAME' = rownames(x),
                                      x)})

        data_IDs <- lapply(X, function(x){merge(x, genes_id)})

        # Add back to MAE object
        metadata(MAE)[["data_IDs"]] <- data_IDs

        return(MAE)

    # For separate analysis
    } else if (method == 's') {

        # For each list within the list add the entrez/ ensembl IDs
        data_IDs <- Map(function(x, y) lapply(
                              x, function(dat) {dat$GENENAME <- row.names(dat);
                              merge(dat, y)}), filtered_genelist, list(miR_IDs,
                                                                       mRNA_IDs
                                                                       ))
        # Add back to MAE object
        metadata(MAE)[["data_IDs"]] <- data_IDs

        return(MAE)

     } else print('Enter c or s as method.')
}
