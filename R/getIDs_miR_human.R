#' @title getIDs_miR_human
#' @description getIDs_miR_human will produce ensembl and entrez data for human
#' microRNAs. It will aslo produce adjusted ensembl and entrez for IDs that are
#'specific to multiple microRNAs.
#' @param MAE MultiAssayObject created by StartObject
#' @param miR Dataframe. Rownames are genes and columns are 
#' results of DE.
#' @return MAE object with 4 more dataframes consisting of ID information.
#' @export
#' @importFrom clusterProfiler bitr
#' @importFrom biomaRt useMart getBM
#' @import org.Hs.eg.db
#' @usage getIDs_miR_human(MAE, miR)
#' @examples
#' library(biomaRt)
#' library(org.Hs.eg.db)
#' library(clusterProfiler)
#' hs_miR -> miR
#' miR[1:100,] -> miR
#' gsub(rownames(miR), pattern = "\\.", replacement =  "-") -> rownames(miR)
#' StartObject(miR = miR, mRNA = NULL) -> MAE
#' getIDs_miR_human(MAE = MAE, miR = MAE@ExperimentList$miR) -> MAE
getIDs_miR_human <- function(MAE, miR){
if(missing(MAE)) stop('Use the StartObject function.');
if(missing(miR)) stop('Add microRNA as.data.frame.');
miR$Genes <- miR$MicroRNA <- rownames(miR)
miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "")
miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "")
miR$MicroRNA <- MicroRNA_full(miRdf = miR$MicroRNA, species = 'hsa')
miR_entrez <- bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', 
toType = 'ENTREZID', OrgDb = org.Hs.eg.db)
miR_ensembl <- bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', 
toType = 'ENSEMBL', OrgDb = org.Hs.eg.db)
miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'MicroRNA',
by.y = 'GENENAME', all = TRUE)
miR_merged <- merge(x = miR_merged, y = miR_entrez, by.x = 'MicroRNA', 
by.y = 'GENENAME', all = TRUE)
miR_merged <- miR_merged[!duplicated(miR_merged$Genes),]
miR_merged <- miR_merged[order(miR_merged$Genes),]
miR_merged$ENTREZID_adjusted <- non_unique(Col = miR_merged$ENTREZID, 
sep = ".", suffix = "")
miR_merged$ENSEMBL_adjusted <- non_unique(Col = miR_merged$ENSEMBL, sep = ".",
suffix = "")
miR_merged <- miR_merged[! duplicated(miR_merged$Genes),]
miR_merged <- miR_merged[order(miR_merged$Genes),]
rownames(miR_merged) <- miR_merged$Genes
MAE@ExperimentList$miR_entrez <- as.data.frame(cbind(GENENAME = 
rownames(miR_merged), ID = miR_merged$ENTREZID))
MAE@ExperimentList$miR_adjusted_entrez <- as.data.frame(cbind(
GENENAME = rownames(miR_merged),ID = miR_merged$ENTREZID_adjusted))
MAE@ExperimentList$miR_ensembl <- as.data.frame(cbind( GENENAME = 
rownames(miR_merged), ID = miR_merged$ENSEMBL))
MAE@ExperimentList$miR_adjusted_ensembl <- as.data.frame(cbind( 
GENENAME = rownames(miR_merged), ID = miR_merged$ENSEMBL_adjusted))
return(MAE)
}
