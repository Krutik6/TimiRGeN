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
#' library(MultiAssayExperiment)
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' miR <- mm_miR
#' miR <- miR[1:100,]
#' MAE <- StartObject(miR = miR, mRNA = NULL)
#' MAE <- getIDs_miR_mouse(MAE, assay(MAE, 1))
getIDs_miR_mouse <- function(MAE, miR){
    if(missing(MAE)) stop('Use the StartObject function.');
    if (missing(miR)) stop('Add microRNA as.data.frame. Rownames are genes and
                            columns are results from differential
                            expression analysis.');
    
    miR$Genes <- miR$names <- rownames(miR)
    miR$Genes <- sub(x = miR$Genes, pattern = "-3p", replacement = "")
    miR$Genes <- sub(x = miR$Genes, pattern = "-5p", replacement = "")
    miR$Genes <- MicroRNA_full(miRdf = miR$Genes, species = 'mmu')
    
    # Get entrez and ensembl IDs
    miR_entrez <- suppressWarnings(bitr(geneID = miR$Genes, 
                                        fromType = 'GENENAME',
                                        toType = 'ENTREZID',
                                        OrgDb = org.Mm.eg.db))
                                   
    miR_ensembl <- suppressWarnings(bitr(geneID = miR$Genes, 
                                         fromType = 'GENENAME',
                                         toType = 'ENSEMBL', 
                                         OrgDb = org.Mm.eg.db))
    
    miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'Genes',
                        by.y = 'GENENAME', all = TRUE)
    
    miR_merged <- merge(x = miR_merged, y = miR_entrez, 
                        by.x = 'Genes',  by.y = 'GENENAME', all = TRUE)
    
    miR_merged <- miR_merged[!duplicated(miR_merged$names),]
    miR_merged <- miR_merged[order(miR_merged$names),]

    # Get adjusted entrezIDs and esembl IDs
    miR_merged$ENTREZID_adjusted <- non_unique(Col = miR_merged$ENTREZID,
                                               sep = ".", suffix = "")
    
    miR_merged$ENSEMBL_adjusted <- non_unique(Col = miR_merged$ENSEMBL, 
                                              sep = ".",  suffix = "")
    
    rownames(miR_merged) <- miR_merged$names
    
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
                                        ID = miR_merged$ENSEMBL_adjusted)))))
    
    MAE <- c(MAE, MAE2)
    
return(MAE)
}
