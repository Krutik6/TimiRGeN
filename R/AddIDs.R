#' @title AddIDs
#' @aliases AddIDs
#' @description Adds entrez or ensembl IDs to the nested dataframes within
#' the filtered_genelist.
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
#' @import BiocManager
#' @usage AddIDs(method, filtered_genelist, miR_IDs, mRNA_IDs)
#' @examples
#'library(clusterProfiler)
#'library(org.Mm.eg.db)
#'mm_miR -> miR
#'mm_mRNA -> mRNA
#'StartObject(miR = miR, mRNA = mRNA) -> MAE
#'getIDs_miR_mouse(MAE = MAE, miR = MAE@ExperimentList$miR) -> MAE
#'getIDs_mRNA_mouse(MAE = MAE, mRNA = MAE@ExperimentList$mRNA) -> MAE
#' 
#'CombineGenes(miR_data = MAE@ExperimentList$miR, 
#'mRNA_data = MAE@ExperimentList$mRNA) -> MAE@ExperimentList$genetic_data
#' 
#'GenesList(method = 'c', genetic_data = MAE@ExperimentList$genetic_data,
#'timeString = 'D') -> MAE@metadata$genelist
#' 
#'SignificantVals(method = "c", geneList = MAE@metadata$genelist, 
#'maxVal = 0.05, stringVal = "adjPVal") -> MAE@metadata$filtered_genelist
#' 
#'AddIDs(method = "c", filtered_genelist = MAE@metadata$filtered_genelist, 
#'miR_IDs = MAE@ExperimentList$miR_entrez,
#'mRNA_IDs = MAE@ExperimentList$mRNA_entrez) -> MAE@metadata$gene_entrez
AddIDs <- function(method, filtered_genelist, miR_IDs, mRNA_IDs){
if (missing(method)) stop('method should be s for separate analysis and
c for combined analysis.')
if (missing(filtered_genelist)) stop('Input list of nested as.data.frames');
if (missing(miR_IDs)) stop('Input miR as.data.frame which contains a list
of genenames and entrezids/ ensembl gene names.');
if (missing(mRNA_IDs)) stop('Input miRNA as.data.frame which contains a
list of genenames and entrezids/ ensembl
gene names.');
if (method == 'c') {
colnames(miR_IDs) <- c("GENENAME", "ID")
colnames(mRNA_IDs) <- c("GENENAME", "ID")
geneIDs <- rbind(miR_IDs, mRNA_IDs)
genes_id <- geneIDs[! duplicated(geneIDs[[1]]),]
X <- lapply(filtered_genelist, function(x){cbind('GENENAME' = rownames(x),
x)})
Y <- lapply(X, function(x){merge(x, genes_id)})
return(Y)
} else if (method == 's') {
colnames(miR_IDs) <- c("GENENAME", "ID")
colnames(mRNA_IDs) <- c("GENENAME", "ID")
X <- Map(function(x, y) lapply(x, function(dat) {dat$GENENAME <-
row.names(dat);
merge(dat, y)}), filtered_genelist, list(miR_IDs, mRNA_IDs))
return(X)
} else {stop('Please insert method c for combined analysis or s
for seperate analysis')}
}
