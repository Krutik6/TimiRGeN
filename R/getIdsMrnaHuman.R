#' @title getIdsMrnaHuman
#' @description getIdsMrnaHuman will produce ensembl and entrez ID dataframes
#' for human mRNAs. These will be stored as 2 individual assays within a MAE.
#' @param MAE MultiAssayExperiment to store the output of getIdsMrnaHuman.
#' It is recommended to use the MAE which contains output from startObject.
#' @param mRNA A Dataframe. Rownames are genes and columns are results of DE.
#' This should be found as an assay within the MAE used in the startObject
#' function. Please read vignette for nomenclature guidance.
#' @param mirror String to identify which biomaRt server is best. This is based
#' on location. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'www'.
#' @return 2 new dataframes. One with entrez information and another with
#' ensembl information. Output will be stored as assays in the input MAE.
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
getIdsMrnaHuman <- function(MAE, mRNA, mirror = "www"){

    if(missing(MAE)) stop('
                          MAE is missing.
                          Add MAE to store output of getIdsMrnaHuman. Please use
                          startObject first.')

    if(missing(mRNA)) stop('
                           mRNA is missing.
                           Add mRNA dataframe. Please use startObject first.
                           The output of startObject will be stored as an
                           assay within the MAE used in the startObject
                           function.')

    mRNA$Genes <- rownames(mRNA)

    # download a human mart
    human <- suppressMessages(biomaRt::useEnsembl(
                                         "ensembl",
                                          dataset="hsapiens_gene_ensembl",
                                          host = paste0(mirror,".ensembl.org")))

    # Get IDs
    glist <- suppressMessages(biomaRt::getBM(
                                      attributes = c("external_gene_name",
                                                     "ensembl_gene_id",
                                                      "entrezgene_id"),
                                      filters = "external_gene_name",
                                      values = mRNA$Genes,
                                      mart = human,
                                      uniqueRows = TRUE))

    # merge acquired data to input data
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
