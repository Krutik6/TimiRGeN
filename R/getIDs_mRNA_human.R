#' @title getIDs_mRNA_human
#' @description getIDs_miR_human will produce ensembl and entrez data
#'for human mRNAs.
#' @param MAE MultiAssayObject created by StartObject
#' @param mRNA Dataframe. Rownames should be genenames and colnames should
#'be results from DE.
#' @param mirror String to identify which biomaRt repo is best for the users
#' locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 2 new dataframes. One with entrez information and another with
#'ensembl information.
#' @importFrom biomaRt useEnsembl
#' @export
#' @usage getIDs_mRNA_human(MAE, mRNA, mirror)
#' @examples
#' library(biomaRt)
#' library(MultiAssayExperiment)
#' mRNA <- mm_mRNA
#' mRNA <- mRNA[1:20,]
#' MAE <- StartObject(miR = NULL, mRNA = mRNA)
#' MAE <- getIDs_mRNA_human(MAE = MAE, mRNA = assay(MAE, 2), mirror = 'useast')
getIDs_mRNA_human <- function(MAE, mRNA, mirror = "useast"){
    
    if(missing(MAE)) stop('Use the StartObject function.');
    if(missing(mRNA)) stop('Add microRNA as.data.frame. Rownames are genes
                            and columns are results from differential
                           expression analysis.');
    
    mRNA$Genes <- rownames(mRNA)
    
    human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
                                 GRCh=37, host = paste0(mirror, ".ensembl.org"))
    
    # Get IDs
    glist <- getBM(attributes = c("external_gene_name", "ensembl_gene_id",
                                  "entrezgene_id"),
                   filters = "external_gene_name", 
                   values = mRNA$Genes,
                   mart = human, uniqueRows = TRUE)
    
    m_dat <- merge(x = mRNA, y = glist, by.x = 'Genes', 
                   by.y = 'external_gene_name', all = TRUE)
    
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
