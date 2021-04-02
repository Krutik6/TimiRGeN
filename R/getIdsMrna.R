#' @title getIdsMrna
#' @description  getIdsMrna will produce ensembl and entrez ID dataframes
#' for mRNAs. These will be stored as 2 individual assays within a MAE.
#' @param MAE MultiAssayExperiment to store the output of getIdsMrna.
#' It is recommended to use the MAE which contains output from startObject.
#' @param mRNA Dataframe. Rownames are genes and columns are results of DE.
#' This should be found as an assay within the MAE used in the startObject
#' function. Please read vignette for nomenclature guidance.
#' @param mirror String to identify which biomaRt server is best. This is based
#' on location. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'www'.
#' @param species Species of interest. E.g. mmusculus or hsapiens.
#' @return 2 new dataframes in the MAE. One with entrez information and
#' the other with ensembl gene ID information.
#' @export
#' @importFrom biomaRt useEnsembl getBM
#' @usage getIdsMrna(MAE, mRNA, mirror, species)
#' @examples
#'data(mm_mRNA)
#'
#' mRNA <- mm_mRNA[1:20,]
#'
#' MAE <- startObject(miR = NULL, mRNA = mRNA)
#'
#' MAE <- getIdsMrna(MAE = MAE, mRNA = assay(MAE, 2), mirror = 'useast',
#'                       species = 'mmusculus')
getIdsMrna <- function(MAE, mRNA, mirror = 'www', species){

    if(missing(MAE)) stop('MAE is missing. Add MAE to store output of getIdsMrna. Please use startObject first.')

    if (missing(mRNA)) stop('mRNA is missing. Add mRNA dataframe. Please use startObject first. The output of startObject will be stored as an assay within the MAE used in the startObject function.')

    if (missing(species)) stop('species is missing. Add species of interest e.g. "mmusculus", "hsapiens"')

    # Get a  mart
    mart <- suppressMessages(biomaRt::useEnsembl("ensembl",
                                                dataset=paste0(species,
                                                               '_gene_ensembl'),
                                                host = paste0(mirror,
                                                              ".ensembl.org")))

    # Get IDs
    glist <- suppressMessages(biomaRt::getBM(attributes = c(
                                                "external_gene_name",
                                                "ensembl_gene_id",
                                                "entrezgene_id"),
                                              filters = "external_gene_name",
                                              values = rownames(mRNA),
                                              mart = mart, uniqueRows = TRUE))

    # Remove duplicated data
    glist <- glist[! duplicated(glist$external_gene_name),]

    gene_data <- cbind(mRNA, rownames(mRNA))


    # Merge retrieved IDs to input data
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
