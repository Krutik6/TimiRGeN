#' @title getIDs_miR_human
#' @description getIDs_miR_human will produce ensembl and entrez data for human
#' microRNAs. It will aslo produce adjusted ensembl and entrez for IDs that are
#'specific to multiple microRNAs.
#' @param miR Dataframe. Rownames are genes and columns are results of DE.
#' @return 4 dataframes consisting of ID information.
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom biomaRt useMart getBM
#' @import org.Hs.eg.db
#' @usage getIDs_miR_human(miR)
#' @examples
#' library(biomaRt)
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' hs_miR -> miR
#' gsub(rownames(miR), pattern = "\\.", replacement =  "-") -> rownames(miR)
#' getIDs_miR_human(miR)
getIDs_miR_human <- function(miR){
if(missing(miR)) stop('Add microRNA dataframe. Rownames are genes and columns
are results from differential expression analysis.')
miR$Genes <- miR$MicroRNA <- rownames(miR)
gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "") -> miR$MicroRNA
gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "") -> miR$MicroRNA
MicroRNA_full(miRdf = miR$MicroRNA, species = 'hsa') -> miR$MicroRNA
bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', toType = 'ENTREZID',
OrgDb = org.Hs.eg.db) -> miR_entrez
bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', toType = 'ENSEMBL',
OrgDb = org.Hs.eg.db) -> miR_ensembl
merge(x = miR, y = miR_ensembl, by.x = 'MicroRNA', by.y = 'GENENAME',
all = TRUE) -> miR_merged
merge(x = miR_merged, y = miR_entrez, by.x = 'MicroRNA', by.y = 'GENENAME',
all = TRUE) -> miR_merged
miR_merged[!duplicated(miR_merged$Genes),] -> miR_merged
miR_merged[order(miR_merged$Genes),] -> miR_merged
non_unique(Col = miR_merged$ENTREZID, sep = ".",
suffix = "") -> miR_merged$ENTREZID_adjusted
non_unique(Col = miR_merged$ENSEMBL, sep = ".",
suffix = "") -> miR_merged$ENSEMBL_adjusted
miR_merged[! duplicated(miR_merged$Genes),] -> miR_merged
miR_merged <- miR_merged[order(miR_merged$Genes),]
rownames(miR_merged) <- miR_merged$Genes
miR_entrez <- as.data.frame(cbind(GENENAME = rownames(miR_merged),
ID = miR_merged$ENTREZID))
miR_adjusted_entrez <- as.data.frame(cbind(GENENAME = rownames(miR_merged),
ID = miR_merged$ENTREZID_adjusted))
miR_ensembl <- as.data.frame(cbind(GENENAME = rownames(miR_merged),
ID = miR_merged$ENSEMBL))
miR_adjusted_ensembl <- as.data.frame(cbind(GENENAME = rownames(miR_merged),
ID = miR_merged$ENSEMBL_adjusted))
miR_nested <- list(miR_entrez = miR_entrez,
miR_adjusted_entrez = miR_adjusted_entrez,
miR_ensembl = miR_ensembl,
miR_adjusted_ensembl = miR_adjusted_ensembl)
return(list2env(miR_nested, .GlobalEnv))
}
