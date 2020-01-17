#' @title getIDs_mRNA_mousetohuman
#' @description getIDs_miR_mousetohuman will produce ensembl and entrez data
#'for mouse mRNAs. Also a new dataframe of human named mRNAs will be made.
#' @param mRNA Dataframe. Rownames should be genenames and colnames should be
#'results from DE.
#' @param mirror String to identify which biomaRt repo is best for the users
#'locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @return 3 dataframes. One with human entrez information, another with human
#'ensembl information and a new dataframe of human named genes with the DE
#'results.
#' @usage getIDs_mRNA_mousetohuman(mRNA, mirror)
#' @export
#' @examples
#' library(biomaRt)
#' mm_mRNA -> mRNA
#' mRNA[1:200,] -> mRNA
#' getIDs_mRNA_mousetohuman(mRNA = mRNA, mirror = 'useast')
getIDs_mRNA_mousetohuman <- function(mRNA, mirror = 'useast'){
        if (missing(mRNA)) stop('Add microRNA dataframe. Rownames are genes and
        columns are results from differential
        expression analysis.')
        mRNA$musGenes <- rownames(mRNA)
        human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
        GRCh=37, host = paste0(mirror, ".ensembl.org"))
        mouse <- biomaRt::useEnsembl("ensembl",
        dataset="mmusculus_gene_ensembl",
        GRCh=37, host = paste0(mirror, ".ensembl.org"))
        genesV2 <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol",
        values = mRNA$musGenes ,
        mart = mouse, attributesL = c("hgnc_symbol",
        "entrezgene_id",
        "ensembl_gene_id"),
        martL = human, uniqueRows=TRUE)
        genesV2[! duplicated(genesV2$MGI.symbol),] -> no_dups_MGI
        no_dups_MGI[! duplicated(no_dups_MGI$HGNC.symbol),] -> no_dups_HGNC
        no_dups_HGNC[order(no_dups_HGNC$MGI.symbol),] -> order_MGI
        order_MGI[which(order_MGI$MGI.symbol %in% rownames(mRNA) == TRUE),
        ] -> Which_MGI_mRNArows
        merge(x = mRNA, y = Which_MGI_mRNArows, by.x = "musGenes",
        by.y = "MGI.symbol") -> Human_merged
        rownames(Human_merged) <- Human_merged$musGenes
        as.data.frame(cbind(GENENAME = rownames(Human_merged),
        ID = Human_merged$EntrezGene.ID)) -> mRNA_entrez
        as.data.frame(cbind(GENENAME =rownames(Human_merged),
        ID = Human_merged$Gene.stable.ID)) -> mRNA_ensembl
        Human_merged$musGenes <- Human_merged$HGNC.symbol <-
        Human_merged$EntrezGene.ID <- Human_merged$Gene.stable.ID <- NULL
        mRNA_nested <- list(mRNA_human_renamed = Human_merged,
        mRNA_entrez = mRNA_entrez,
        mRNA_ensembl = mRNA_ensembl)
return(list2env(mRNA_nested, .GlobalEnv))
}
