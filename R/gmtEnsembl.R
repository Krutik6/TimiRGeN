#' @title gmtEnsembl
#' @description Change entrez IDs in path_gene and path_data into ensembl IDs.
#' @param MAE MultiAssayExperiment where the ensembl wikipathways data can
#' be stored.
#' @param path_gene path_gene data from downloadGMT.
#' @param path_data path_data from downloadGMT.
#' @param orgDB organism Library e.g. org.Hs.eg.db
#' @return 2 dataframes. Wikipathway data based on ensembl IDs.
#' @export
#' @importFrom clusterProfiler bitr
#' @usage gmtEnsembl(MAE, path_gene, path_data, orgDB)
#' @examples
#' library(org.Mm.eg.db)
#'
#' miR <- mm_miR
#'
#' mRNA <- mm_mRNA
#'
#' MAE <- startObject(miR = miR, mRNA = mRNA)
#'
#' MAE <- dloadGmt(MAE, speciesInitial = "Mm")
#'
#' MAE <- gmtEnsembl(MAE = MAE, assay(MAE, 3),
#'                    assay(MAE, 5), org.Mm.eg.db)
gmtEnsembl <- function(MAE, path_gene, path_data, orgDB){

    if(missing(MAE)) stop('Add an MAE object.')

    if(missing(path_gene)) stop('Add gene data from dloadGmt')

    if(missing(path_data)) stop('Add path data from dloadGmt')

    if(missing(orgDB)) stop('Add org.xx.eg.db, either Hs or Mm')

    path_gene <- as.data.frame(path_gene)

    path_data <- as.data.frame(path_data)

    # Retreive ensembl IDs
    ent_ens <- suppressWarnings(clusterProfiler::bitr(geneID = path_gene$gene,
                                                      fromType = 'ENTREZID',
                                                      toType = 'ENSEMBL',
                                                      OrgDb = orgDB))

    # Add ensembl IDs to path_gen and path_data
    path_gene_ensembl <- merge(x = path_gene, y = ent_ens, by.x = 'gene',
                               by.y = 'ENTREZID')

    path_data_ensembl <- merge(x = path_gene_ensembl,
                               y = path_data,
                               by.x = 'gene',
                               by.y = 'gene')

    # Edit path_data
    path_data_ensembl_corrected <- path_data_ensembl[path_data_ensembl[,
                                                  2] == path_data_ensembl[,4], ]

    # Remove uneeded columns
    path_gene_ensembl$gene <- NULL

    path_data_ensembl_corrected$gene <- NULL

    path_data_ensembl_corrected$wpid.y <- NULL

    # Make sure second column in each is 'gene'
    names(path_gene_ensembl)[2] <- 'gene'

    names(path_data_ensembl_corrected)[2] <- 'gene'

    # Save to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                'path_gene_ensembl' = path_gene_ensembl,
                'path_data_ensembl' = path_data_ensembl_corrected
    )))

    MAE <- c(MAE, MAE2)
return(MAE)
}
