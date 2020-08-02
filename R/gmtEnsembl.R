#' @title gmtEnsembl
#' @description Change entrez IDs in path_gene and path_data into ensembl IDs.
#' Will create two new dataframes with ensembl IDs and wikipathway information,
#' and store these in as assays in a MAE object.
#' @param MAE MultiAssayExperiment where the ensembl wikipathways data can
#' be stored. It is recommended to use the same MAE object which contains
#' output from dloadGmt.
#' @param path_gene Dataframe with wpid and entrezgene IDs.
#' path_gene is from the dloadGmt function. It will be stored as an assay
#' in the MAE used in the dloadGmt function.
#' @param path_data Dataframe with wpid, wikipathway names and entrezgene IDs.
#' path_data is from the dloadGmt function. It will be stored as an assay
#' in the MAE used in the dloadGmt function.
#' @param orgDB Load an appropriate org Library e.g. org.Hs.eg.db.
#' @return 2 dataframes. Wikipathway data and associated gene names based on
#' ensembl IDs.
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

    if(missing(MAE)) stop('Add a MAE object to store output from gmtEnsebl.
                          Please use dloadGmt first.')

    if(missing(path_gene)) stop('Add dataframe with entrezgene IDs and
                                wikipathways IDs. Please use the dloadGmt
                                function first. The output of dloadGmt should
                                be stored as an assay within the MAE used in
                                the dloatGmt function.')

    if(missing(path_data)) stop('Add dataframe with entrezgene IDs,
                                wikipathway IDs and wikipathways names.
                                Please use the dloadGmt function first.
                                The output of dloadGmt should be stored as an
                                assay within the MAE used in the dloatGmt
                                function.')

    if(missing(orgDB)) stop('Please load either org.Hs.eg.db or org.Mm.eg.db.')

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

    # Remove unneeded columns
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
