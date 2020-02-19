#' @title getIDs_miR_mouse
#' @description getIDs_miR_human will produce ensembl and entrez data for mouse
#'microRNAs. It will aslo produce adjusted ensembl and entrez for IDs that
#'are specific to multiple microRNAs.
#' @param MAE MultiAssayObject created by StartObject
#' @param miR Dataframe. Rownames are genes and columns are results of DE.
#' @return MAE object with 4 more dataframes consisting of ID information.
#' @import org.Mm.eg.db
#' @export
#' @usage getIDs_miR_mouse(MAE, miR)
#' @examples
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' mm_miR -> miR
#' miR[1:100,] -> miR
#' StartObject(miR = miR, mRNA = NULL) -> MAE
#' getIDs_miR_mouse(MAE, MAE@ExperimentList$miR) -> MAE
getIDs_miR_mouse <- function(MAE, miR){
if (missing(miR)) stop('Add microRNA as.data.frame. Rownames are genes and
columns are results from differential expression analysis.')
miR$Genes <- miR$names <- rownames(miR)
miR$Genes <- sub(x = miR$Genes, pattern = "-3p", replacement = "")
miR$Genes <- sub(x = miR$Genes, pattern = "-5p", replacement = "")
miR$Genes <- MicroRNA_full(miRdf = miR$Genes, species = 'mmu')
miR_entrez <- bitr(geneID = miR$Genes, fromType = 'GENENAME',
toType = 'ENTREZID', OrgDb = org.Mm.eg.db)
miR_ensembl <- bitr(geneID = miR$Genes, fromType = 'GENENAME',
toType = 'ENSEMBL', OrgDb = org.Mm.eg.db)
miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'Genes', by.y = 'GENENAME',
all = TRUE)
miR_merged <- merge(x = miR_merged, y = miR_entrez, by.x = 'Genes', 
by.y = 'GENENAME', all = TRUE)
miR_merged <- miR_merged[!duplicated(miR_merged$names),]
miR_merged <- miR_merged[order(miR_merged$names),]
miR_merged$ENTREZID_adjusted <- non_unique(Col = miR_merged$ENTREZID, sep = ".",
suffix = "")
miR_merged$ENSEMBL_adjusted <- non_unique(Col = miR_merged$ENSEMBL, sep = ".",
suffix = "")
rownames(miR_merged) <- miR_merged$names
MAE@ExperimentList$miR_entrez <- as.data.frame(cbind(GENENAME = 
rownames(miR_merged),ID = miR_merged$ENTREZID))
MAE@ExperimentList$miR_adjusted_entrez <- as.data.frame(
cbind(GENENAME = rownames(miR_merged),
ID = miR_merged$ENTREZID_adjusted))
MAE@ExperimentList$miR_ensembl <- as.data.frame(cbind(
GENENAME = rownames(miR_merged),
ID = miR_merged$ENSEMBL))
MAE@ExperimentList$miR_adjusted_ensembl <- as.data.frame(cbind(
GENENAME = rownames(miR_merged),
ID = miR_merged$ENSEMBL_adjusted))
return(MAE)
}
