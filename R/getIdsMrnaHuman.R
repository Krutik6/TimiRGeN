#' @title getIdsMrnaHuman
#' @description getIDs_miR_human will produce ensembl and entrez data
#'for human mRNAs.
#' @param MAE MultiAssayObject created by StartObject
#' @param mRNA Dataframe. Rownames should be genenames and colnames should
#'be results from DE and time points.
#' @param mirror String to identify which biomaRt repo is best for the users
#' locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 2 new dataframes. One with entrez information and another with
#'ensembl information.
#' @export
#' @importFrom biomaRt useEnsembl getBM
#' @usage getIdsMrnaHuman(MAE, mRNA, mirror)
#' @examples
#' mRNA <- mm_mRNA
#'
#' mRNA <- mRNA[1:20,]
#'
#' MAE <- startObject(miR = NULL, mRNA = mRNA)
#'
#' MAE <- getIdsMrnaHuman(MAE = MAE, mRNA = assay(MAE, 2), mirror = 'useast')
getIdsMrnaHuman <- function(MAE, mRNA, mirror = "useast"){

    if(missing(MAE)) stop('Use the startObject function.')

    if(missing(mRNA)) stop('Add mRNA data frame.')

    mRNA$Genes <- rownames(mRNA)

    # download a human mart
    human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
                                 host = paste0(mirror, ".ensembl.org"))

    # Get IDs
    glist <- biomaRt::getBM(attributes = c("external_gene_name",
                                           "ensembl_gene_id",
                                            "entrezgene_id"),
                            filters = "external_gene_name",
                            values = mRNA$Genes,
                            mart = human, uniqueRows = TRUE)

    # merge aquired data to input data
    m_dat <- merge(x = mRNA, y = glist, by.x = 'Genes',
                   by.y = 'external_gene_name', all = TRUE)

    # remove duplicates and order data
    m_dat <- m_dat[! duplicated(m_dat$Genes),]

    m_dat <- m_dat[order(m_dat$Genes),]

    rownames(m_dat) <- m_dat$Genes

    # Save to MAE object
    MAE2 <- suppressMessages(MultiAssayExperiment(list(
                                                mRNA_entrez = data.frame(cbind(
                                                    GENENAME = rownames(m_dat),
                                                    ID = m_dat$entrezgene_id)),
                                                mRNA_enembl = data.frame(cbind(
                                                    GENENAME = rownames(m_dat),
                                                    ID = m_dat$ensembl_gene_id)
                                                    ))))

    MAE <- c(MAE, MAE2)
return(MAE)
}
