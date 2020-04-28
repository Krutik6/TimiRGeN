#' @title enrichWiki
#' @aliases enrichWiki
#' @description Finds which wikipathways are enriched within the data.
#' @param MAE MultiAssayExperiment object.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param ID_list List of ensembl or entrez IDs for each sample.
#' @param orgDB Library of species specific data.
#' @param path_gene genes-wikipathwayIDs dataframe. From downladGMT or
#' GMT_ensembl.
#' @param path_name wikipathwayIDs-wikipathway fullnames dataframe.
#' From downladGMT.
#' @param ID Either ENTREZID or ENSEMBL.
#' @param universe A column of gene IDs to be used as the background for
#' gene set enrichment. To find which pathways are most enriched use
#' path_gene$gene
#' @param pvalcutoff Default is 0.05. P value cutuff point.
#' @param qvaluecutoff Default is 0.2. q value cutoff point.
#' @param padjustmethod Default is 'BH'. Which pvalue adjustment method should
#' be used. Look into clusterProfiler function enricher for more options.
#' @return A large list of data which identifies which wikipathways are
#' most enriched in the input data which can be stored in an MAE object.
#' @export
#' @importFrom clusterProfiler enricher setReadable
#' @usage enrichWiki(MAE, method = '', ID_list, orgDB, path_gene, path_name,
#'                   ID = '', universe, pvalcutoff, qvaluecutoff,
#'                   padjustmethod)
#' @examples
#' library(org.Mm.eg.db)
#' miR <- mm_miR
#' mRNA <- mm_mRNA
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' metadata(MAE)[["e_list"]] <- e_list
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- enrichWiki(MAE = MAE, method = 'c', ID_list = metadata(MAE)[[1]],
#'                   orgDB = org.Mm.eg.db, path_gene = assay(MAE, 3),
#'                   path_name = assay(MAE, 4), ID = "ENTREZID",
#'                   universe = assay(MAE, 3)[[2]])
enrichWiki <- function(MAE, method, ID_list, orgDB, path_gene, path_name, ID,
                       universe = universe, pvalcutoff = 0.05,
                       qvaluecutoff = 0.2, padjustmethod = "BH"){

    if (missing(MAE)) stop('Check if you are using the correct MAE object.');
    if (missing(ID_list)) stop('Input list of entrezIDs.');
    if (missing(orgDB)) stop('Input an org.Db package.');
    if (missing(path_gene)) stop ('Use downloadGMT to get wikipathway data.');
    if (missing(path_name)) stop ('Use downloadGMT to get wikipathway data.');
    if (missing(ID)) stop ('Enter either ENTREZID or ENSEMBL. Should match
                            wpid data used.');

    if (missing(universe)) stop ('Add a list of gene IDs to be used as the
                                  background for gene set enrichment.
                                  To find which pathways are most enriched use
                                  path_gene$gene / assay(MAE, i)[[2]]');

    metadata <- `metadata<-` <- NULL

    # Perform enricher function on each list in the list
    lst2 <- lapply(ID_list, function(x){
    enricher(x, TERM2GENE = path_gene,
             TERM2NAME = path_name,
             universe = universe,
             pvalueCutoff = pvalcutoff,
             qvalueCutoff = qvaluecutoff,
             pAdjustMethod = padjustmethod)})

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
seperate analysis.')}
}
