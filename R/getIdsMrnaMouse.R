#' @title getIdsMrnaMouse
#' @description  getIDs_miR_mouse will produce ensembl and entrez data for
#' mouse mRNAs.
#' @param MAE MultiAssayObject created by StartObject
#' @param mRNA Dataframe. Rownames should be genenames and colnames should be
#'results from DE.
#' @param mirror tring to identify which biomaRt repo is best for the users
#' locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 2 new dataframes in the MAE object. One with entrez information and
#' the other with ensembl information.
#' @export
#' @usage getIdsMrnaMouse(MAE, mRNA, mirror)
#' @examples
#' mRNA <- mm_mRNA
#' mRNA <- mRNA[1:20,]
#' MAE <- startObject(miR = NULL, mRNA = mRNA)
#' MAE <- getIdsMrnaMouse(MAE = MAE, mRNA = assay(MAE, 2), mirror = 'useast')
getIdsMrnaMouse <- function(MAE, mRNA, mirror = 'useast'){

    if(missing(MAE)) stop('Use the StartObject function.');
    if (missing(mRNA)) stop('Add microRNA as.data.frame. Rownames are genes and
                            columns are results from differential
                            expression analysis.');

    mouse <- biomaRt::useEnsembl("ensembl",dataset="mmusculus_gene_ensembl",
                                 GRCh=37, host = paste0(mirror, ".ensembl.org"))

    # Get IDs
    glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                                  "entrezgene_id"),
                                  filters = "external_gene_name",
                                  values = rownames(mRNA),
                                  mart = mouse, uniqueRows = TRUE)

    glist <- glist[! duplicated(glist$external_gene_name),]

    gene_data <- cbind(mRNA, rownames(mRNA))

    m_dat <- merge(x = gene_data, y = glist, by.x = 'rownames(mRNA)',
                   by.y = 'external_gene_name', all = TRUE)

    rownames(m_dat) <- m_dat$`rownames(mRNA)`

    # Save to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                            mRNA_entrez = data.frame(cbind(
                                                GENENAME = rownames(m_dat),
                                                ID = m_dat$entrezgene_id)),
                                            mRNA_enembl = data.frame(cbind(
                                                GENENAME = rownames(m_dat),
                                                ID = m_dat$ensembl_gene_id)))))

    MAE <- c(MAE, MAE2)

return(MAE)
}
