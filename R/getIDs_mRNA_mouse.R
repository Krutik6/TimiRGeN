#' @title getIDs_mRNA_mouse
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
#' @usage getIDs_mRNA_mouse(MAE, mRNA, mirror)
#' @examples
#' library(biomaRt)
#' mm_mRNA -> mRNA
#' mRNA[1:20,] -> mRNA
#' StartObject(miR = NULL, mRNA = mRNA) -> MAE
#' getIDs_mRNA_mouse(MAE = MAE, mRNA = MAE@ExperimentList$mRNA, 
#' mirror = 'useast') -> MAE
getIDs_mRNA_mouse <- function(MAE, mRNA, mirror = 'useast'){
if (missing(mRNA)) stop('Add microRNA as.data.frame. Rownames are genes and
columns are results from differential
expression analysis.')
mouse <- biomaRt::useEnsembl("ensembl",dataset="mmusculus_gene_ensembl",
GRCh=37, host = paste0(mirror, ".ensembl.org"))
glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
"entrezgene_id"),
filters = "external_gene_name", values = rownames(mRNA),
mart = mouse, uniqueRows = TRUE)
glist <- glist[! duplicated(glist$external_gene_name),]
gene_data <- cbind(mRNA, rownames(mRNA))
m_dat <- merge(x = gene_data, y = glist, by.x = 'rownames(mRNA)', by.y =
'external_gene_name', all = TRUE)
rownames(m_dat) <- m_dat$`rownames(mRNA)`
MAE@ExperimentList$mRNA_entrez <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$entrezgene_id)) 
MAE@ExperimentList$mRNA_ensembl <- as.data.frame(cbind(GENENAME = rownames(
m_dat), ID = m_dat$ensembl_gene_id))
return(MAE)
}