#' @title gmtEnsembl
#' @description Change entrez IDs in path_gene and path_data into ensembl IDs.
#' Will create two new dataframes with ensembl IDs and wikipathway information.
#' @param MAE MultiAssayExperiment which will store the output of gmtEnsembl.
#' It is recommended to use the same MAE object which contains output from
#' dloadGmt.
#' @param path_gene Dataframe with wikipathway IDs and entrezgene IDs.
#' path_gene is from the dloadGmt function. It will be stored as an assay
#' within the MAE used in the dloadGmt function.
#' @param path_data Dataframe with wikipathway IDs, wikipathway names and
#' entrezgene IDs. path_data is from the dloadGmt function. It will be stored
#' as an assay within the MAE used in the dloadGmt function.
#' @param orgDB Load the appropriate db package e.g. org.Hs.eg.db if human
#' wikipathways are being used.
#' @return 2 dataframes. One containing wikipathway IDs and ensembl gene IDs,
#' and the other containing  wikipathway IDs, ensembl gene IDs and wikipathway
#' names. Output will be stored as assays in the input MAE.
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

    if(missing(MAE)) stop('
                          MAE is missing.
                          Add a MAE to store output from gmtEnsembl.
                          Please use dloadGmt first.')

    if(missing(path_gene)) stop('
                                path_gene is missing.
                                Add dataframe with entrezgene IDs and
                                wikipathways IDs. Please use the dloadGmt
                                function first. The output of dloadGmt should
                                be stored as an assay within the MAE used in
                                the dloadGmt function.')

    if(missing(path_data)) stop('
                                path_data is missing.
                                Add dataframe with entrezgene IDs,
                                wikipathway IDs and wikipathways names.
                                Please use the dloadGmt function first.
                                The output of dloadGmt should be stored as an
                                assay within the MAE used in the dloadGmt
                                function.')

    if(missing(orgDB)) stop('
                             orgDB is missing.
                             Please load either org.Hs.eg.db or org.Mm.eg.db.')

    path_gene <- as.data.frame(path_gene)

    path_data <- as.data.frame(path_data)

    # Retreive ensembl IDs
    ent_ens <- suppressMessages(suppressWarnings(
                                 clusterProfiler::bitr(geneID = path_gene$gene,
                                 fromType = 'ENTREZID',
                                 toType = 'ENSEMBL',
                                 OrgDb = orgDB)))

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
