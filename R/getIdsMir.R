#' @title getIdsMir
#' @description getIdsMir will produce ensembl and entrez ID data for
#' microRNAs. It will also produce adjusted ensembl and entrez for IDs
#' that are specific to microRNAs that share an ID. They will be stored as
#' 4 individual assays in a MAE. org.Mm.eg.db must be loaded prior to using
#' this function.
#' @param MAE MultiAssayExperiment to store the output of getIdsMir.
#' It is recommended to use the MAE which contains output from startObject.
#' @param miR A Dataframe. Rownames are genes and columns are results of DE.
#' This should be found as an assay within the MAE used in the startObject
#' function. Please read vignette for nomenclature guidance.
#' @param orgDB org.xx.eg.db package which corresponds to the species being
#' analysed.
#' @param miRPrefix microRNA prefix for the species being analysed e.g. 'mmu',
#' 'hsa', 'rno' ect.
#' @return 4 dataframes consisting of either entrez or ensembl ID information.
#' 2 of these will be adjusted for shared IDs. Output will be stored as assays
#' in the input MAE.
#' @export
#' @importFrom clusterProfiler bitr
#' @usage getIdsMir(MAE, miR, orgDB, miRPrefix)
#' @examples
#' library(org.Mm.eg.db)
#'
#' data(mm_miR)
#'
#' # Make sure miRNA gene name nomenclature is correct for TimiRGeN analysis!
#'
#' miR <- mm_miR[1:100,]
#'
#' MAE <- startObject(miR = miR, mRNA = NULL)
#'
#' MAE <- getIdsMir(MAE, assay(MAE, 1), orgDB = org.Mm.eg.db, miRPrefix = 'mmu')
getIdsMir <- function(MAE, miR, orgDB, miRPrefix){

    if(missing(MAE)) stop('MAE is missing. Add MAE to store output of getIdsMir. Please use startObject first.');

    if (missing(miR)) stop('miR is missing. Add microRNA dataframe. Please use startObject first. The output of startObject will be stored as an assay within the MAE used in the startObject function.')

    if(missing(orgDB)) stop('orgDB is missing. Add org.xx.eg.db package which is relevant for the analysis.');

    if(missing(miRPrefix)) stop('miRPrefix is missing. Add microRNA prefix e.g. "mmu" for mouse or "rno" for rat.');

        miR$Genes <- miR$MicroRNA <- rownames(miR)

        # remove -3p and -5p for now
        miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "")

        miR$MicroRNA <- gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "")

        # Use micrornaFull to standardise miRNA names
        miR$MicroRNA <- micrornaFull(miRdf = miR$MicroRNA, species = miRPrefix)

        # Get entrez and ensembl IDs
        miR_entrez <- suppressMessages(suppressWarnings(clusterProfiler::bitr(
                                                        geneID = miR$MicroRNA,
                                                        fromType = 'GENENAME',
                                                        toType = 'ENTREZID',
                                                        OrgDb = orgDB)))

        miR_ensembl <- suppressMessages(suppressWarnings(clusterProfiler::bitr(
                                                        geneID = miR$MicroRNA,
                                                        fromType = 'GENENAME',
                                                        toType = 'ENSEMBL',
                                                        OrgDb = orgDB)))
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
