#' @title getIdsMirHuman
#' @description getIdsMirHuman will produce ensembl and entrez ID data for
#' human microRNAs. It will also produce adjusted ensembl and entrez for IDs
#' that are specific to microRNAs that share an ID. They will be stored as
#' 4 individual assays in a MAE. org.Hs.eg.db must be loaded prior to use of
#' this function.
#' @param MAE MultiAssayExperiment to store the output of getIdsMirHuman.
#' It is recommended to use the MAE which contains output from startObject.
#' @param miR Dataframe. Rownames are genes and columns are results of DE.
#' This should be found as an assay within the MAE used in the
#' startObject function. Please read vignette for nomenclature guidance.
#' @return MAE object with 4 more dataframes consisting of ID information.
#' @export
#' @importFrom clusterProfiler bitr
#' @usage getIdsMirHuman(MAE, miR)
#' @examples
#' library(org.Hs.eg.db)
#'
#' miR <- hs_miR
#'
#' # Make sure miRNA gene name nomenclature is correct for TimiRGeN analysis!
#'
#' rownames(miR) <- gsub(rownames(miR), pattern = "\\.", replacement = "-")
#'
#' rownames(miR) <- sub("-$", "*", rownames(miR))
#'
#' miR <- miR[1:100,]
#'
#' MAE <- startObject(miR = miR, mRNA = NULL)
#'
#' MAE <- getIdsMirHuman(MAE, assay(MAE, 1))
getIdsMirHuman <- function(MAE, miR){

    if(missing(MAE)) stop('MAE object to store output of getIdsMirHuman.
                          Please use startObject first.')

    if(missing(miR)) stop('Add microRNA dataframe. Please use startObject
                           first. The output of startObject will be stored as an
                           assay within the MAE used in the startObject
                           function.')

    miR$Genes <- miR$MicroRNA <- rownames(miR)

    # remove -3p and -5p for now
    miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "")

    miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "")

    # Use micrornaFull to standardise miRNA names
    miR$MicroRNA <- micrornaFull(miRdf = miR$MicroRNA, species = 'hsa')

    # Get entrez and ensembl IDs
    miR_entrez <- suppressMessages(suppressWarnings(clusterProfiler::bitr(
                                                        geneID = miR$MicroRNA,
                                                        fromType = 'GENENAME',
                                                        toType = 'ENTREZID',
                                                        OrgDb = org.Hs.eg.db)))

    miR_ensembl <- suppressMessages(suppressWarnings(clusterProfiler::bitr(
                                                         geneID = miR$MicroRNA,
                                                         fromType = 'GENENAME',
                                                         toType = 'ENSEMBL',
                                                         OrgDb = org.Hs.eg.db)))

    # Merge ensembl and entrez data
    miR_merged <- merge(x = miR, y = miR_ensembl, by.x = 'MicroRNA',
                        by.y = 'GENENAME', all = TRUE)

    miR_merged <- merge(x = miR_merged, y = miR_entrez, by.x = 'MicroRNA',
                        by.y = 'GENENAME', all = TRUE)

    # remove duplices and order
    miR_merged <- miR_merged[!duplicated(miR_merged$Genes),]

    miR_merged <- miR_merged[order(miR_merged$Genes),]

    # Get adjusted entrez and ensembl IDs using nonUnique
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
