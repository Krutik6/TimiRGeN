#' @title getIdsMirHuman
#' @description getIdsMirHuman will produce ensembl and entrez data for human
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
#' @usage getIdsMirHuman(MAE, miR)
#' @examples
#' library(org.Hs.eg.db)
#' miR <- hs_miR
#' rownames(miR) <- gsub(rownames(miR), pattern = "\\.", replacement = "-")
#' rownames(miR) <- sub("-$", "*", rownames(miR))
#' miR <- miR[1:100,]
#' MAE <- startObject(miR = miR, mRNA = NULL)
#' MAE <- getIdsMirHuman(MAE, assay(MAE, 1))
getIdsMirHuman <- function(MAE, miR){
    if(missing(MAE)) stop('Use the StartObject function.');
    if(missing(miR)) stop('Add microRNA as.data.frame.');

    miR$Genes <- miR$MicroRNA <- rownames(miR)
    miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "")
    miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "")
    miR$MicroRNA <- micrornaFull(miRdf = miR$MicroRNA, species = 'hsa')

    # Get entrez and ensembl IDs
    miR_entrez <- suppressWarnings(bitr(geneID = miR$MicroRNA,
                                        fromType = 'GENENAME',
                                        toType = 'ENTREZID',
                                        OrgDb = org.Hs.eg.db))

    miR_ensembl <- suppressWarnings(bitr(geneID = miR$MicroRNA,
                                         fromType = 'GENENAME',
                                         toType = 'ENSEMBL',
                                         OrgDb = org.Hs.eg.db))

    miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'MicroRNA',
                        by.y = 'GENENAME', all = TRUE)

    miR_merged <- merge(x = miR_merged, y = miR_entrez, by.x = 'MicroRNA',
                        by.y = 'GENENAME', all = TRUE)

    miR_merged <- miR_merged[!duplicated(miR_merged$Genes),]
    miR_merged <- miR_merged[order(miR_merged$Genes),]
    # Get adjusted entrez and ensembl IDs
    miR_merged$ENTREZID_adjusted <- nonUnique(Col = miR_merged$ENTREZID,
                                               sep = ".", suffix = "")

    miR_merged$ENSEMBL_adjusted <- nonUnique(Col = miR_merged$ENSEMBL,
                                              sep = ".",
                                              suffix = "")

    miR_merged <- miR_merged[! duplicated(miR_merged$Genes),]
    miR_merged <- miR_merged[order(miR_merged$Genes),]

    rownames(miR_merged) <- miR_merged$Genes
    # Save to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                        miR_entrez = data.frame(cbind(
                                            GENENAME = rownames(miR_merged),
                                            ID = miR_merged$ENTREZID)),
                                        miR_ensembl = data.frame(cbind(
                                            GENENAME = rownames(miR_merged),
                                            ID = miR_merged$ENSEMBL)),
                                        miR_adjusted_entrez = data.frame(cbind(
                                            GENENAME = rownames(miR_merged),
                                            ID = miR_merged$ENTREZID_adjusted)),
                                        miR_adjusted_ensembl = data.frame(cbind(
                                            GENENAME = rownames(miR_merged),
                                            ID = miR_merged$ENSEMBL_adjusted))
                                        )))
    MAE <- c(MAE, MAE2)

return(MAE)
}
