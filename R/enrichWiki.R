#' @title enrichWiki
#' @description Finds which wikipathways are enriched within the data. This
#' function uses gene set enrichment analysis from clusterProfiler to
#' find enriched signalling pathways. Each time point is analysed
#' individually. In the case of separated TimiRGeN analysis, each gene type and
#' time point are analysed individually.
#' @param MAE MultiAssayExperiment which will store the output from enrichWiki.
#' It is recommended to use the MAE object which stores the output from the
#' dloadGmt function.
#' @param method Either 'c' or 's', respectively for combined or separated
#' analysis.
#' @param ID_list List of ensembl or entrez IDs for each sample.
#' This is the output from eNames function. This will be found as
#' metadata within the MAE used in the eNames function.
#' @param orgDB DB package of the species being analysed. e.g. org.Mm.eg.db
#' if mouse miR-mRNA data is being looked into.
#' @param path_gene Dataframe containing pathway ID - gene ID information. This
#' is output from either dloadGmt or gmtEnsembl. It will be stored as an assay
#' within the MAE used in dloadGmt or gmtEnsembl.
#' @param path_name Dataframe containing pathway ID - pathway names information.
#' This is output from dloadGmt. It will be stored as an assay within the MAE
#' used in dloadGmt.
#' @param ID Either "ENTREZID" or "ENSEMBL". This should be the same as the ID
#' type used for ID_list. dloadGmt loads data as entrez gene IDs and gmtEnsembl
#' converts this to ensembl gene IDs.
#' @param universe A column of gene IDs to be used as the background for
#' gene set enrichment. IDs should be stored as characters. It is recommended to
#' use all genes found within the wikipathways of the species being analysed as
#' background i.e. path_gene$gene or universe = assay(MAE, i)[[2]]/ MAE[[i]][2].
#' To add a unique universe, create a list of gene IDs (entrezID or ensembl)
#' which are classed as characters.
#' @param pvalcutoff Default is 0.05. P value cut-off point.
#' @param qvaluecutoff Default is 0.2. q value cut-off point.
#' @param padjustmethod Default is 'BH'. This sets the pvalue adjustment
#' method. Look into the enricher function from clusterProfiler for more info.
#' @return A large list which identifies which wikipathways are most enriched
#' at each time point of the input data. Output will be stored as metadata in
#' the input MAE.
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
#' metadata(MAE)[["e_list"]] <- e_list_mouse
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

    if (missing(MAE)) stop('
                           MAE is missing.
                           Add MAE. Results from enrichWiki will be
                           stored in this MAE. Please use eNames, dloadGmt /
                           gmtEnsembl first.')

    if (missing(ID_list)) stop('
                               ID_list is missing.
                               Add list of entrezIDs. Please use the eNames
                               function fist. The output of this should be in
                               the metadata the MAE used in the eNames function.
                               ')

    if (missing(orgDB)) stop('
                             orgDB is missing.
                             Please load org.Mm.eg.db or org.Hs.eg.db before
                             running this function.')

    if (missing(path_gene)) stop ('
                                  path_gene is missing.
                                  Please use dloadGmt/ gmtEnsembl first to get a
                                  wikipathway ID - gene ID dataframe. Output of
                                  these functions are stored as assays within
                                  the MAE used in dloadGmt/ gmtEnsembl functions
                                  .')

    if (missing(path_name)) stop ('
                                  path_name is missing.
                                  Please use dloadGmt first to get wikipathway
                                  ID - pathway names dataframe. Will be stored
                                  as an assay in the MAE object used in the
                                  dloadGmt function.')

    if (missing(ID)) stop ('
                           ID is missing.
                           Add either "ENTREZID" or "ENSEMBL". Should match ID
                           types found in ID_list and path_gene parameters.')

    if (missing(universe)) stop ('
                                  universe is missing.
                                  Add a list of gene IDs to be used as the
                                  background for gene set enrichment.
                                  The recommended method is to use
                                  path_gene$gene / assay(MAE, i)[[2]].
                                  For alternative analysis, add own list of gene
                                  IDs.')

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

    } else print('Enter c or s as method.')
}
