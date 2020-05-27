#' @title getIdsMirMouse
#' @description getIDs_miR_human will produce ensembl and entrez data for mouse
#'microRNAs. It will aslo produce adjusted ensembl and entrez for IDs that
#'are specific to microRNAs that share an ID.
#' @param MAE MultiAssayObject created by StartObject
#' @param miR Dataframe. Rownames are genes and columns are results of DE.
#' @return MAE object with 4 more dataframes consisting of ID information.
#' @export
#' @importFrom clusterProfiler bitr
#' @usage getIdsMirMouse(MAE, miR)
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' miR <- miR[1:100,]
#'
#' MAE <- startObject(miR = miR, mRNA = NULL)
#'
#' MAE <- getIdsMirMouse(MAE, assay(MAE, 1))
getIdsMirMouse <- function(MAE, miR){

    if(missing(MAE)) stop('Use the startObject function.')

    if (missing(miR)) stop('Add microRNA as.data.frame. Rownames are genes and
                            columns are results from differential
                            expression analysis.')

    miR$Genes <- miR$names <- rownames(miR)

    # remove -5p or -3p data
    miR$Genes <- sub(x = miR$Genes, pattern = "-3p", replacement = "")

    miR$Genes <- sub(x = miR$Genes, pattern = "-5p", replacement = "")

    # use micrornaFull to standardise names
    miR$Genes <- micrornaFull(miRdf = miR$Genes, species = 'mmu')

    # Get entrez and ensembl IDs
    miR_entrez <- suppressWarnings(clusterProfiler::bitr(
                                        geneID = miR$Genes,
                                        fromType = 'GENENAME',
                                        toType = 'ENTREZID',
                                        OrgDb = org.Mm.eg.db))

    miR_ensembl <- suppressWarnings(clusterProfiler::bitr(
                                         geneID = miR$Genes,
                                         fromType = 'GENENAME',
                                         toType = 'ENSEMBL',
                                         OrgDb = org.Mm.eg.db))

    # merge retreived data
    miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'Genes',
                        by.y = 'GENENAME', all = TRUE)

    miR_merged <- merge(x = miR_merged, y = miR_entrez,
                        by.x = 'Genes',  by.y = 'GENENAME', all = TRUE)

    # remove duplicated data and order data
    miR_merged <- miR_merged[!duplicated(miR_merged$names),]

    miR_merged <- miR_merged[order(miR_merged$names),]

    # Get adjusted entrezIDs and esembl IDs
    miR_merged$ENTREZID_adjusted <- nonUnique(Col = miR_merged$ENTREZID,
                                               sep = ".", suffix = "")

    miR_merged$ENSEMBL_adjusted <- nonUnique(Col = miR_merged$ENSEMBL,
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
