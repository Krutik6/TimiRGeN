#' @title enrichWiki
#' @description Finds which wikipathways are enriched within the data. This
#' function uses gene set enrichment analysis from clusterProfiler to
#' find enriched pathways from wikipathways. Each time point is analysed
#' individually. In the case of TimiRGeN analysis using the 's' analysis,
#' each gene type and time point is analysed individually. Data is stored as
#' metadata of a MAE.
#' @param MAE MultiAssayExperiment. Results of enrichWiki will be
#' stored in the metadata. It is recommended to use the MAE object created
#' by dloadGmt.
#' @param method Either 'c' or 's' for combined or separated analysis.
#' @param ID_list List of ensembl or entrez IDs for each sample.
#' This is the output from eNames function. This will be found within the
#' metadata of the MAE used in the eNames function.
#' @param orgDB Library of species specific data. e.g. org.Mm.eg.db.
#' @param path_gene path_gene dataframe. From dloadGmt or gmtEnsembl
#' function, and will be stored as assays within the MAE used in dloadGmt or
#' gmtEnsembl.
#' @param path_name path_name dataframe from dloadGmt. Will be stored as an
#' assay within the MAE used in dloadGmt.
#' @param ID Either "ENTREZID" or "ENSEMBL". This should be the same as the ID
#' type used for ID_list. dloadGMT loads data as ENTREZID, so use gmtEnsembl
#' function to get ENSEMBL data.
#' @param universe A column of gene IDs to be used as the background for
#' gene set enrichment. Recommended use is to use all genes found within the
#' wikipathways of the species of interest as background i.e.
#' path_gene$gene or universe = assay(MAE, i)[[2]]. To add unique universe,
#' create a list of gene IDs (entrezID or ensembl) to contrast against.
#' @param pvalcutoff Default is 0.05. P value cutuff point.
#' @param qvaluecutoff Default is 0.2. q value cutoff point.
#' @param padjustmethod Default is 'BH'. This sets the  pvalue adjustment
#' method to be used. Look into clusterProfiler function enricher for more
#' options.
#' @return A large list of data which identifies which wikipathways are
#' most enriched in the input data which can be stored in an MAE object.
#' @export
#' @importFrom clusterProfiler enricher setReadable
#' @usage enrichWiki(MAE, method = '', ID_list, orgDB, path_gene, path_name,
#'                   ID = '', universe, pvalcutoff, qvaluecutoff,
#'                   padjustmethod)
#' @examples
#' library(org.Mm.eg.db)
#'
#' MAE <- MultiAssayExperiment()
#'
#' metadata(MAE)[["e_list"]] <- e_list
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- enrichWiki(MAE = MAE, method = 'c', ID_list = metadata(MAE)[[1]],
#'                   orgDB = org.Mm.eg.db, path_gene = assay(MAE, 1),
#'                   path_name = assay(MAE, 2), ID = "ENTREZID",
#'                   universe = assay(MAE, 1)[[2]])
enrichWiki <- function(MAE, method, ID_list, orgDB, path_gene, path_name, ID,
                       universe = universe, pvalcutoff = 0.05,
                       qvaluecutoff = 0.2, padjustmethod = "BH"){

    if (missing(MAE)) stop('Add MAE object. Results from enrichWiki will be
                           stored in this MAE. Please use eNames, dloadGmt /
                           gmtEnsembl first.')

    if (missing(ID_list)) stop('Add list of entrezIDs. Please use
                               eNames function fist. The output of this
                               should be in the metadata the MAE used in the
                               eNames function.')

    if (missing(orgDB)) stop('Please load org.Mm.eg.db or org.Hs.eg.db before
                             running this function.')

    if (missing(path_gene)) stop ('Please use dloadGmt/ gmtEnsembl to get
                                  wikipathway data. Output of these functions
                                  are stored as assays within MAE objects used
                                  in dloadGmt/ gmtEnsembl functions.')

    if (missing(path_name)) stop ('Please use dloadGMT to get wikipathway data.
                                  Will be stored as an assay in the MAE object
                                  used in the dloadGmt function.')

    if (missing(ID)) stop ('Enter either ENTREZID or ENSEMBL. Should match ID
                           types used in eNames output and wpid data.')

    if (missing(universe)) stop ('Add a list of gene IDs to be used as the
                                  background for gene set enrichment.
                                  The recommended method is to use
                                  path_gene$gene / assay(MAE, i)[[2]].
                                  For alternative analysis, a user can
                                  add their own list of gene IDs.')

    metadata <- `metadata<-` <- NULL

    # Perform enricher function on each list in the list
    lst2 <- lapply(ID_list, function(x){
                                    enricher(x, TERM2GENE = path_gene,
                                             TERM2NAME = path_name,
                                             universe = universe,
                                             pvalueCutoff = pvalcutoff,
                                             qvalueCutoff = qvaluecutoff,
                                             pAdjustMethod = padjustmethod)})

    # Extract information in an organised  way
    sigwiki <- lapply(lst2, function(x){setReadable(x, orgDB, keyType = ID)})

    # If c look into first level of names() to get new names
    if (method == 'c') {

        names(sigwiki) <- paste0(names(sigwiki), '_wikipathways', sep='')

        metadata(MAE)[["sigwiki"]] <- sigwiki

    return(MAE)

    # If s look into second level of names() to get new names
    }else if (method == 's') {

        names(sigwiki) <- gsub(x = names(sigwiki), pattern = '\\.',
                               replacement = '_wikipathways')

        metadata(MAE)[["sigwiki"]] <- sigwiki

    return(MAE)

}else {stop('Please insert method c for combined analysis or s for
            separate analysis.')}
}
