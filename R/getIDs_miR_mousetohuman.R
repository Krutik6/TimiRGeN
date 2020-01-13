#' @title getIDs_miR_mousetohuman
#' @description getIDs_miR_mousetohuman will produce ensembl and entrez data
#' for mouse microRNAs. It will aslo produce adjusted ensembl and entrez for
#' IDs that are specific to multiple microRNAs. Finally a new miR_human_named
#' file will be made.
#' @param miR Dataframe. Rownames should be genenames and colnames should be
#' results from DE.
#' @param mirror String to identify which biomaRt repo is best for the users
#'locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @importFrom biomaRt getLDS
#' @return 5 dataframes consisting of ID information, and a new dataframe of
#' human named miRs with the DE results.
#' @export
#' @examples
#' library(biomaRt)
#' library(clusterProfiler)
#' library(org.Hs.eg.db)
#' mm_miR -> miR
#' miR[1:100,] -> miR
#' getIDs_miR_mousetohuman(miR = miR, mirror = 'useast')
getIDs_miR_mousetohuman <- function(miR, mirror = 'useast'){
if (missing(miR)) stop('Add mouse microRNA dataframe. Rownames are genes and
columns are results from differential expression analysis.')
cbind(miR, Gene = rownames(miR), name = rownames(miR)) -> miR
gsub(x = miR$Gene, pattern = "-3p", replacement = "") -> miR$Gene
gsub(x = miR$Gene, pattern = "-5p", replacement = "") -> miR$Gene
human = useEnsembl(biomart = "ensembl",
dataset = "hsapiens_gene_ensembl", mirror = mirror)
mouse = useEnsembl(biomart = "ensembl",
dataset = "mmusculus_gene_ensembl", mirror = mirror)
bmt = getLDS(attributes = c("mirbase_id"), filters = "mirbase_id",
values = miR$Gene, mart = mouse, attributesL = c("mirbase_id"), martL = human)
colnames(bmt) <- c("Mm_n", "Hs_n")
gsub(bmt$Mm_n, pattern = "mir", replacement = "miR") -> bmt$Mm_n
gsub(bmt$Hs_n, pattern = "mir", replacement = "miR") -> bmt$Hs_n
merge(x = miR, y = bmt, by.x = "Gene", by.y = "Mm_n", all = TRUE) -> mh_m
mh_m[!duplicated(mh_m$name),] -> mh_d
rownames(mh_d) <- mh_d$name
mh_d$Hs_n[is.na(mh_d$Hs_n)] <- as.character(mh_d$Gene[is.na(mh_d$Hs_n)])
gsub(mh_d$Hs_n, pattern = 'mmu', replacement = 'hsa') -> mh_d$Hs_n
mh_d[order(mh_d$name),] -> mh_o
MicroRNA_full(mh_o$Hs_n, 'hsa') -> mh_o$microRNA
bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENTREZID",
OrgDb = org.Hs.eg.db) -> Y
merge(x =mh_o, y=Y, by.x="microRNA", by.y="GENENAME", all=TRUE)->miR_en
bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENSEMBL",
OrgDb = org.Hs.eg.db) -> Z
merge(x=miR_en, y=Z, by.x="microRNA", by.y="GENENAME", all=TRUE) -> miR_IDs
non_unique(Col = miR_IDs$Hs_n, sep = '-', suffix = 'p') -> miR_IDs$Hs_n
gsub(miR_IDs$Hs_n, pattern = "-1p", replacement = "-3p") -> miR_IDs$Hs_n
gsub(miR_IDs$Hs_n, pattern = "-2p", replacement = "-5p") -> miR_IDs$Hs_n
miR_IDs$ENTREZID -> miR_IDs$ENTREZID_adj
non_unique(Col=miR_IDs$ENTREZID_adj,sep ='.',suffix='')->miR_IDs$ENTREZID_adj
miR_IDs$ENSEMBL -> miR_IDs$ENSEMBL_adj
non_unique(Col=miR_IDs$ENSEMBL_adj,sep ='.',suffix='')->miR_IDs$ENSEMBL_adj
miR_IDs[! duplicated(miR_IDs$Hs_n),] -> miR_IDs
miR_IDs$Hs_n -> rownames(miR_IDs)
as.data.frame(cbind(GENENAME=rownames(miR_IDs),ID=miR_IDs$ENTREZID))->miR_entrez
as.data.frame(cbind(GENENAME = rownames(miR_IDs),
ID = miR_IDs$ENTREZID_adj)) -> miR_adjusted_entrez
as.data.frame(cbind(GENENAME=rownames(miR_IDs),ID=miR_IDs$ENSEMBL))->miR_ensembl
as.data.frame(cbind(GENENAME = rownames(miR_IDs),
ID = miR_IDs$ENSEMBL_adj)) -> miR_adjusted_ensembl
miR_IDs$microRNA<-miR_IDs$Gene<-miR_IDs$name<-miR_IDs$Hs_n<-miR_IDs$ENTREZID<-
miR_IDs$ENSEMBL<-miR_IDs$ENTREZID_adj<-miR_IDs$ENSEMBL_adj<-NULL
miR_nested <- list(miR_human_renamed = miR_IDs, miR_entrez = miR_entrez,
miR_adjusted_entrez = miR_adjusted_entrez, miR_ensembl = miR_ensembl,
miR_adjusted_ensembl = miR_adjusted_ensembl)
return(list2env(miR_nested, .GlobalEnv))
}
