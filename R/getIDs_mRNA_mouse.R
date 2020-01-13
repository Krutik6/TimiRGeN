#' @title getIDs_mRNA_mouse
#' @description  getIDs_miR_mouse will produce ensembl and entrez data for
#' mouse mRNAs.
#' @param mRNA Dataframe. Rownames should be genenames and colnames should be
#'results from DE.
#' @param mirror tring to identify which biomaRt repo is best for the users
#' locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 2 new dataframes. One with entrez information and another with
#' ensembl information.
#' @export
#' @usage getIDs_mRNA_mouse(mRNA, mirror)
#' @examples
#' library(biomaRt)
#' mm_mRNA -> mRNA
#' mRNA[1:200,] -> mRNA
#' getIDs_mRNA_mouse(mRNA)
getIDs_mRNA_mouse <- function(mRNA, mirror = 'useast'){
        if (missing(mRNA)) stop('Add microRNA dataframe. Rownames are genes and
        columns are results from differential
        expression analysis.')
        mouse = useEnsembl(biomart = "ensembl",
        dataset = "mmusculus_gene_ensembl", mirror = mirror)
        glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
        "entrezgene_id"),
        filters = "external_gene_name", values = rownames(mRNA),
        mart = mouse, uniqueRows = TRUE)
        glist[! duplicated(glist$external_gene_name),] -> glist
        cbind(mRNA, rownames(mRNA)) -> gene_data
        merge(x = gene_data, y = glist, by.x = 'rownames(mRNA)', by.y =
        'external_gene_name', all = TRUE) -> m_dat
        rownames(m_dat) <- m_dat$`rownames(mRNA)`
        as.data.frame(cbind(GENENAME = rownames(m_dat),
        ID = m_dat$entrezgene_id)) -> mRNA_entrez
        as.data.frame(cbind(GENENAME = rownames(m_dat),
        ID = m_dat$ensembl_gene_id)) -> mRNA_ensembl
        mRNA_nested <- list(mRNA_entrez = mRNA_entrez,
        mRNA_ensembl = mRNA_ensembl)
return(list2env(mRNA_nested, .GlobalEnv))
}
