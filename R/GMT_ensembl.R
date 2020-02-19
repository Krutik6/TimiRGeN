#' @title GMT_ensembl
#' @description Change entrez IDs in wpid2gene and wp_data into ensembl IDs
#' @param MAE MultiAssayExperiment where the ensembl wikipathways data can 
#' be stored.
#' @param path_gene wikipathway-gene data from downloadGMT.
#' @param path_data wikipathway-gene-wikipathway full name data from
#' downloadGMT.
#' @param orgDB organism Library e.g. org.Hs.eg.db
#' @return 2 dataframes. Wikipathway data based on ensembl IDs.
#' @export
#' @usage GMT_ensembl(MAE, path_gene, path_data, orgDB)
#' @examples
#' library(org.Mm.eg.db)
#' mm_miR -> miR
#' mm_mRNA -> mRNA
#' StartObject(miR = miR, mRNA = mRNA) -> MAE
#' 
#' dloadGMT(MAE = MAE, speciesInitials = "Mm") -> MAE
#' 
#' GMT_ensembl(MAE = MAE, MAE@ExperimentList$path_gene,
#' MAE@ExperimentList$path_data, org.Mm.eg.db) -> MAE
GMT_ensembl <- function(MAE, path_gene, path_data, orgDB){
path_gene <- as.data.frame(path_gene)
path_data <- as.data.frame(path_data)
ent_ens <- bitr(geneID = path_gene$gene, fromType = 'ENTREZID', 
toType = 'ENSEMBL',OrgDb = orgDB)
path_gene_ensembl <- merge(x = path_gene, y = ent_ens, by.x = 'gene',
by.y = 'ENTREZID') 
path_data_ensembl <- merge(x = path_gene_ensembl, y = path_data, by.x = 'gene',
by.y = 'gene') 
path_data_ensembl_corrected <- path_data_ensembl[path_data_ensembl[,
2] == path_data_ensembl[,4], ]
path_gene_ensembl$gene <- NULL
path_data_ensembl_corrected$gene <- NULL
path_data_ensembl_corrected$wpid.y <- NULL
names(path_gene_ensembl)[2] <- 'gene'
names(path_data_ensembl_corrected)[2] <- 'gene'
MAE@ExperimentList$path_data_ensembl <- as.data.frame(
path_data_ensembl_corrected)
MAE@ExperimentList$path_gene_ensembl <- as.data.frame(path_gene_ensembl)
return(MAE)
}
