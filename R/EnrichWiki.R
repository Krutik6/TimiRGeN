#' @title EnrichWiki
#' @aliases EnrichWiki
#' @description Finds which wikipathways are enriched within the data.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param e_list List of ensembl or entrez IDs for each sample.
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
#' @usage EnrichWiki(method = '', e_list, orgDB, path_gene, path_name, ID = '',
#' universe, pvalcutoff, qvaluecutoff, padjustmethod)
#' @examples
#' library(org.Mm.eg.db)
#' library(clusterProfiler)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' 
#' e_list -> MAE@metadata$elist
#' dloadGMT(MAE, speciesInitial = "Mm") -> MAE
#' 
#' MAE@metadata$sigwiki <- EnrichWiki(method = "c",
#' e_list = MAE@metadata$elist,
#' orgDB = org.Mm.eg.db, 
#' path_gene = MAE@ExperimentList$path_gene, 
#' path_name = MAE@ExperimentList$path_name, 
#' ID = "ENTREZID", 
#' universe = MAE@ExperimentList$path_gene$gene)
#' 
EnrichWiki <- function(method, e_list, orgDB, path_gene, path_name, ID,
universe, pvalcutoff = 0.05, qvaluecutoff = 0.2,
padjustmethod = "BH"){
if (missing(e_list)) stop('Input list of entrezIDs.');
if (missing(orgDB)) stop('Input an org.Db package.');
if (missing(path_gene)) stop ('Use downloadGMT to get wikipathway data.');
if (missing(path_name)) stop ('Use downloadGMT to get wikipathway data.');
if (missing(ID)) stop ('Enter either ENTREZID or ENSEMBL. Should match
wpid data used.');
if (missing(universe)) stop ('Add a list of gene IDs to be used as the
background for gene set enrichment.
To find which pathways are most enriched use
path_gene$gene');
lapply(e_list, function(x){
enricher(x, TERM2GENE = path_gene, TERM2NAME = path_name,
universe = universe,
pvalueCutoff = pvalcutoff, qvalueCutoff = qvaluecutoff,
pAdjustMethod = padjustmethod)}) -> lst2
lapply(lst2, function(x){
setReadable(x, orgDB, keyType = ID)}) -> W_list
if (method == 'c') {
paste0(names(W_list), '_wikipathways', sep='') -> names(W_list)
return(W_list)
}else if (method == 's') {
gsub(x = names(W_list), pattern = '\\.',
replacement = '_wikipathways') -> names(W_list)
return(W_list)
}else {stop('Please insert method c for combined analysis or s for
seperate analysis.')}
}