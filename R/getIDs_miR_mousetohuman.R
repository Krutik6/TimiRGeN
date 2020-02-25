#' @title getIDs_miR_mousetohuman
#' @description getIDs_miR_mousetohuman will produce ensembl and entrez data
#' for mouse microRNAs. It will aslo produce adjusted ensembl and entrez for
#' IDs that are specific to multiple microRNAs. Finally a new miR_human_named
#' file will be made.
#' @param MAE MultiAssayObject created by StartObject
#' @param miR Dataframe. Rownames should be genenames and colnames should be
#' results from DE.
#' @param mirror String to identify which biomaRt repo is best for the users
#'locations. Either 'useast', 'uswest', 'asia' or 'www'. Default is 'useast'.
#' @importFrom biomaRt getLDS
#' @return 5 dataframes consisting of ID information, and a new dataframe of
#' human named miRs with the DE results.
#' @export
#' @usage getIDs_miR_mousetohuman(MAE, miR, mirror)
#' @examples
#'library(biomaRt)
#'library(clusterProfiler)
#'library(org.Hs.eg.db)
#'mm_miR -> miR
#'StartObject(miR = miR, mRNA = NULL) -> MAE
#'miR[1:10,] -> miR
#'getIDs_miR_mousetohuman(MAE, miR = MAE@ExperimentList$miR, 
#'mirror = 'useast') -> MAE
getIDs_miR_mousetohuman <- function(MAE, miR, mirror = 'useast'){
miR <-cbind(miR, Gene = rownames(miR), name = rownames(miR))
miR$Gene <- gsub(x = miR$Gene, pattern = "-3p", replacement = "") 
miR$Gene <- gsub(x = miR$Gene, pattern = "-5p", replacement = "")
human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
GRCh=37, host = paste0(mirror, ".ensembl.org"))
mouse <- biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl",
GRCh=37, host = paste0(mirror, ".ensembl.org"))
bmt <- getLDS(attributes = c("mirbase_id"), filters = "mirbase_id",
values = miR$Gene,mart = mouse,attributesL = c("mirbase_id"), martL = human)
colnames(bmt) <- c("Mm_n", "Hs_n")
bmt$Mm_n <- gsub(bmt$Mm_n, pattern = "mir", replacement = "miR")
bmt$Hs_n <- gsub(bmt$Hs_n, pattern = "mir", replacement = "miR")
mh_m <- merge(x = miR, y = bmt, by.x = "Gene", by.y = "Mm_n", all = TRUE) 
mh_d <- mh_m[!duplicated(mh_m$name),] 
rownames(mh_d) <- mh_d$name
mh_d$Hs_n[is.na(mh_d$Hs_n)] <- as.character(mh_d$Gene[is.na(mh_d$Hs_n)])
mh_d$Hs_n <- gsub(mh_d$Hs_n, pattern = 'mmu', replacement = 'hsa')
mh_o <- mh_d[order(mh_d$name),]
mh_o$microRNA <- MicroRNA_full(mh_o$Hs_n, 'hsa')
Y <- bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENTREZID",
OrgDb = org.Hs.eg.db)
Y2 <- Y[! duplicated(Y$GENENAME),] 
miR_en <- merge(x =mh_o, y=Y2, by.x="microRNA", by.y="GENENAME", all=TRUE)
Z <- bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENSEMBL",
OrgDb = org.Hs.eg.db)
Z2 <- Z[! duplicated(Z$GENENAME),]
miR_IDs <- merge(x=miR_en, y=Z2, by.x="microRNA", by.y="GENENAME", all=TRUE) 
miR_IDs$Hs_n <- non_unique(Col = miR_IDs$Hs_n, sep = '-', suffix = 'p')
miR_IDs$Hs_n <- gsub(miR_IDs$Hs_n, pattern = "-1p", replacement = "-3p")
miR_IDs$Hs_n <- gsub(miR_IDs$Hs_n, pattern = "-2p", replacement = "-5p")
miR_IDs$ENTREZID_adj <- miR_IDs$ENTREZID
miR_IDs$ENTREZID_adj <- non_unique(Col=miR_IDs$ENTREZID_adj,sep ='.',suffix='')
miR_IDs$ENSEMBL_adj <- miR_IDs$ENSEMBL
miR_IDs$ENSEMBL_adj <- non_unique(Col=miR_IDs$ENSEMBL_adj,sep ='.',suffix='')
miR_IDs <- miR_IDs[! duplicated(miR_IDs$Hs_n),]
rownames(miR_IDs) <- miR_IDs$Hs_n
MAE@ExperimentList$miR_entrez <- as.data.frame(cbind(GENENAME=rownames(
miR_IDs),ID=miR_IDs$ENTREZID))
MAE@ExperimentList$miR_ensembl <- as.data.frame(cbind(GENENAME=rownames(
miR_IDs),ID=miR_IDs$ENSEMBL))
MAE@ExperimentList$miR_adjusted_entrez <- as.data.frame(cbind(
GENENAME = rownames(miR_IDs), ID = miR_IDs$ENTREZID_adj))
MAE@ExperimentList$miR_adjusted_ensembl <- as.data.frame(cbind(
GENENAME = rownames(miR_IDs), ID = miR_IDs$ENSEMBL_adj))
miR_IDs$microRNA<-miR_IDs$Gene<-miR_IDs$name<-miR_IDs$Hs_n<-miR_IDs$ENTREZID<-
miR_IDs$ENSEMBL<-miR_IDs$ENTREZID_adj<-miR_IDs$ENSEMBL_adj<-NULL
MAE@ExperimentList$miR_human_renamed <- miR_IDs
return(MAE)
}
