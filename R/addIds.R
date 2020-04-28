#' @title addIds
#' @aliases addIds
#' @description Adds entrez or ensembl IDs to the nested dataframes within
#' the filtered_genelist.
#' @param MAE MultiAssayExperiment object.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param filtered_genelist A list of nested dataframes if 'c' or A list of
#' lists with nested dataframes if 's'.
#' @param miR_IDs miR_ensembl or miR_entrez. Use getIDs function to acquire
#' this.
#' @param mRNA_IDs mRNA_ensembl or mRNA_entrez. Use getIDs function to acquire
#' this.
#' @return list of dataframes with entrezIDs and genenames additional
#' as columns which can be stored in metadata section of an MAE.
#' @export
#' @usage addIds(MAE, method, filtered_genelist, miR_IDs, mRNA_IDs)
#' @examples
#'library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#' mRNA <- mm_mRNA
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
addIds <- function(MAE, method, filtered_genelist, miR_IDs, mRNA_IDs){

    if (missing(MAE)) stop('Check if you are using the correct MAE object.');
    if (missing(method)) stop('method should be s for separate analysis and
                                c for combined analysis.');
    if (missing(filtered_genelist)) stop('Input list of nested as.data.frames');
    if (missing(miR_IDs)) stop('Input miR as.data.frame which contains a list
                                of genenames and entrezids/ ensembl gene names.'
                                );
    if (missing(mRNA_IDs)) stop('Input miRNA as.data.frame which contains a
                                list of genenames and entrezids/ ensembl gene
                                names.');

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
        merge(dat, y)}), filtered_genelist, list(miR_IDs, mRNA_IDs))
        # Add back to MAE object
        metadata(MAE)[["data_IDs"]] <- data_IDs
        return(MAE)

    } else {stop('Please insert method c for combined analysis or s
    for seperate analysis')}
}
