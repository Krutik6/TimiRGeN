#' @title getIDs_miR_mouse
#' @description getIDs_miR_human will produce ensembl and entrez data for mouse
#'microRNAs. It will aslo produce adjusted ensembl and entrez for IDs that
#'are specific to multiple microRNAs.
#' @param miR Dataframe. Rownames are genes and columns are results of DE.
#' @return 4 dataframes consisting of ID information.
#' @import org.Mm.eg.db
#' @export
#' @usage getIDs_miR_mouse(miR)
#' @examples
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' mm_miR -> miR
#' getIDs_miR_mouse(miR)
getIDs_miR_mouse <- function(miR){
if (missing(miR)) stop('Add microRNA dataframe. Rownames are genes and
columns are results from differential expression analysis.')
miR$Genes <- miR$names <- rownames(miR)
sub(x = miR$Genes, pattern = "-3p", replacement = "") -> miR$Genes
sub(x = miR$Genes, pattern = "-5p", replacement = "") -> miR$Genes
MicroRNA_full(miRdf = miR$Genes, species = 'mmu') -> miR$Genes
bitr(geneID = miR$Genes, fromType = 'GENENAME', toType = 'ENTREZID',
OrgDb = org.Mm.eg.db) -> miR_entrez
bitr(geneID = miR$Genes, fromType = 'GENENAME', toType = 'ENSEMBL',
OrgDb = org.Mm.eg.db) -> miR_ensembl
merge(x = miR, y = miR_ensembl, by.x = 'Genes', by.y = 'GENENAME',
all = TRUE) -> miR_merged
merge(x = miR_merged, y = miR_entrez, by.x = 'Genes', by.y = 'GENENAME',
all = TRUE) -> miR_merged
miR_merged[!duplicated(miR_merged$names),] -> miR_merged
miR_merged[order(miR_merged$names),] -> miR_merged
non_unique(Col = miR_merged$ENTREZID, sep = ".",
suffix = "") -> miR_merged$ENTREZID_adjusted
non_unique(Col = miR_merged$ENSEMBL, sep = ".",
suffix = "") -> miR_merged$ENSEMBL_adjusted
miR_entrez <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENTREZID))
miR_adjusted_entrez <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENTREZID_adjusted))
miR_ensembl <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENSEMBL))
miR_adjusted_ensembl <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENSEMBL_adjusted))
miR_nested <- list(miR_entrez = miR_entrez,
miR_adjusted_entrez = miR_adjusted_entrez,
miR_ensembl = miR_ensembl,
miR_adjusted_ensembl = miR_adjusted_ensembl)
return(list2env(miR_nested, .GlobalEnv))
}
