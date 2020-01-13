#' @title GMT_ensembl
#' @description Change entrez IDs in wpid2gene and wp_data into ensembl IDs
#' @param path_gene wikipathway-gene data from downloadGMT.
#' @param path_data wikipathway-gene-wikipathway full name data from
#' downloadGMT.
#' @param orgDB organism Library e.g. org.Hs.eg.db
#' @return 2 dataframes. Wikipathway data based on ensembl IDs.
#' @export
#' @usage GMT_ensembl(path_gene, path_data, orgDB)
#' @examples
#' library(org.Mm.eg.db)
#' downloadGMT(speciesInitial = "Mm")
#' GMT_ensembl(path_gene, path_data, org.Mm.eg.db)
#' file.remove("mus.gmt")
GMT_ensembl <- function(path_gene = path_gene, path_data = path_data, orgDB){
        bitr(geneID = path_gene$gene, fromType = 'ENTREZID', toType = 'ENSEMBL',
        OrgDb = orgDB) -> ent_ens
        merge(x = path_gene, y = ent_ens, by.x = 'gene',
        by.y = 'ENTREZID') -> path_gene_ensembl
        merge(x = path_gene_ensembl, y = path_data, by.x = 'gene',
        by.y = 'gene') -> path_data_ensembl
        path_data_ensembl[path_data_ensembl[,2] == path_data_ensembl[
        ,4], ] -> path_data_ensembl_corrected
        path_gene_ensembl$gene <- NULL
        path_data_ensembl_corrected$gene <-
        path_data_ensembl_corrected$wpid.y <- NULL
        names(path_gene_ensembl)[2] <- 'gene'
        names(path_data_ensembl_corrected)[2] <- 'gene'
        list(path_gene_ensembl = path_gene_ensembl,
        path_data_ensembl = path_data_ensembl_corrected) -> path_list_ensembl
return(list2env(path_list_ensembl, .GlobalEnv))
}
