#' @title getIDs_mRNA_human
#' @description getIDs_miR_human will produce ensembl and entrez data
#'for human mRNAs.
#' @param mRNA Dataframe. Rownames should be genenames and colnames should
#'be results from DE.
#' @param mirror String to identify which biomaRt repo is best for the users
#' locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 2 new dataframes. One with entrez information and another with
#'ensembl information.
#' @importFrom biomaRt useEnsembl
#' @export
#' @usage getIDs_mRNA_human(mRNA, mirror)
#' @examples
#' library(biomaRt)
#' hs_mRNA -> mRNA
#' mRNA[1:200,] -> mRNA
#' getIDs_mRNA_human(mRNA, mirror = "useast")
getIDs_mRNA_human <- function(mRNA, mirror = "useast"){
        if(missing(mRNA)) stop('Add microRNA dataframe. Rownames are genes
        and columns are results from differential
        expression analysis.')
        rownames(mRNA) -> mRNA$Genes
        human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
        GRCh=37, host = paste0(mirror, ".ensembl.org"))
        glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
        "entrezgene_id"),
        filters = "external_gene_name", values = mRNA$Genes,
        mart = human, uniqueRows = TRUE)
        merge(x = mRNA, y = glist, by.x = 'Genes',
        by.y = 'external_gene_name', all = TRUE) -> m_dat
        m_dat[! duplicated(m_dat$Genes),] -> m_dat
        m_dat[order(m_dat$Genes),] -> m_dat
        rownames(m_dat) <- m_dat$Genes
        as.data.frame(cbind(GENENAME = rownames(m_dat),
        ID = m_dat$entrezgene_id)) -> mRNA_entrez
        as.data.frame(cbind(GENENAME = rownames(m_dat),
        ID = m_dat$ensembl_gene_id)) -> mRNA_ensembl
        mRNA_nested <- list(mRNA_entrez = mRNA_entrez,
        mRNA_ensembl = mRNA_ensembl)
return(list2env(mRNA_nested, .GlobalEnv))
}
