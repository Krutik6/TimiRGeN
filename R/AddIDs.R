#' @title AddIDs
#' @aliases AddIDs
#' @description Adds entrez or ensembl IDs to the nested dataframes within
#' the filtered_genelist.
#' @param method Respectively either 'c' or 's' for combined or separated
#' analysis.
#' @param filtered_genelist A list of nested dataframes if 'c' or A list of
#' lists with nested dataframes if 's'.
#' @param miR_IDs miR_ensembl or miR_entrez. Use getIDs function to acquire
#' this.
#' @param mRNA_IDs mRNA_ensembl or mRNA_entrez. Use getIDs function to acquire
#' this.
#' @return list of dataframes with entrezIDs and genenames additional
#' as columns.
#' @export
#' @import BiocManager
#' @usage AddIDs(method, filtered_genelist, miR_IDs, mRNA_IDs)
#' @examples
#' library(clusterProfiler)
#' library(org.Mm.eg.db)
#' miR <- mm_miR
#' miR <- miR[1:100,]
#' mRNA <- mm_mRNA
#' mRNA <- mRNA[1:200,]
#' getIDs_miR_mouse(miR = miR)
#' getIDs_mRNA_mouse(mRNA = mRNA)
#' CombineGenes(miR_data = miR, mRNA_data = mRNA) -> genetic_data
#' GenesList(method = 'c', genetic_data = genetic_data,
#' timeString = 'D') -> genelist
#' SignificantVals(method = 'c', geneList = genelist, maxVal = 0.05,
#' stringVal = 'adjPVal') -> filtered_genelist
#' AddIDs(method = 'c', filtered_genelist = filtered_genelist,
#' miR_IDs = miR_entrez, mRNA_IDs = mRNA_entrez) -> gene_entrez
AddIDs <- function(method, filtered_genelist, miR_IDs, mRNA_IDs){
        if (missing(method)) stop('method should be s for separate analysis and
        c for combined analysis.')
        if (missing(filtered_genelist)) stop('Input list of nested dataframes');
        if (missing(miR_IDs)) stop('Input miR dataframe which contains a list
        of genenames and entrezids/ ensembl gene names.');
        if (missing(mRNA_IDs)) stop('Input miRNA dataframe which contains a
        list of genenames and entrezids/ ensembl
        gene names.');
        if (method == 'c') {
        colnames(miR_IDs) <- c("GENENAME", "ID")
        colnames(mRNA_IDs) <- c("GENENAME", "ID")
        rbind(miR_IDs, mRNA_IDs) -> geneIDs
        geneIDs[! duplicated(geneIDs$GENENAME),] -> genes_id
        lapply(filtered_genelist, function(x){cbind('GENENAME' = rownames(x),
        x)})-> X
        lapply(X, function(x){merge(x, genes_id)}) -> Y
        return(Y)
        } else if (method == 's') {
        colnames(miR_IDs) <- c("GENENAME", "ID")
        colnames(mRNA_IDs) <- c("GENENAME", "ID")
        Map(function(x, y) lapply(x, function(dat) {dat$GENENAME <-
        row.names(dat);
        merge(dat, y)}), filtered_genelist, list(miR_IDs, mRNA_IDs)) -> X
        return(X)
        } else {stop('Please insert method c for combined analysis or s
        for seperate analysis')}
}
