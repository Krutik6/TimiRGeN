#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(biomaRt)
library(clusterProfiler)
library(org.Hs.eg.db)
library(MultiAssayExperiment)
#load data
mm_miR -> miR
miR[1:20,]-> miR
MAE <- StartObject(miR = miR, mRNA = NULL)
getIDs_miR_mousetohuman(MAE, MAE@ExperimentList$miR) -> MAE
#check 1
test_that("output of getIDs_miR_mousetohuman is as expected", {
expect_equal(length(rownames(MAE@ExperimentList$miR)),
length(rownames(MAE@ExperimentList$miR_ensembl)))
expect_equal(length(rownames(MAE@ExperimentList$miR_ensembl)),
length(rownames(MAE@ExperimentList$miR_human_renamed)))
expect_equal(MAE@ExperimentList$miR_ensembl$GENENAME,
MAE@ExperimentList$miR_adjusted_ensembl$GENENAME)
expect_equal(MAE@ExperimentList$miR_entrez$GENENAME,
MAE@ExperimentList$miR_adjusted_entrez$GENENAME)
expect_false(isTRUE(all.equal(MAE@ExperimentList$miR_ensembl$ID, 
MAE@ExperimentList$miR_adjusted_ensembl$ID)))
expect_false(isTRUE(all.equal(MAE@ExperimentList$miR_entrez$ID,
MAE@ExperimentList$miR_adjusted_entrez$ID)))
})
#internal checks
miR <-cbind(miR, Gene = rownames(miR), name = rownames(miR))
miR$Gene <- gsub(x = miR$Gene, pattern = "-3p", replacement = "") 
miR$Gene <- gsub(x = miR$Gene, pattern = "-5p", replacement = "")
#check 2
test_that("musGenes is different from rownames", {
expect_equal(length(colnames(miR)), 12)
})
#continue
human <- biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl",
GRCh=37, host = paste0(mirror, ".ensembl.org"))
mouse <- biomaRt::useEnsembl("ensembl", dataset="mmusculus_gene_ensembl",
GRCh=37, host = paste0(mirror, ".ensembl.org"))
bmt <- getLDS(attributes = c("mirbase_id"), filters = "mirbase_id",
values = miR$Gene,mart = mouse,attributesL = c("mirbase_id"), martL = human)
colnames(bmt) <- c("Mm_n", "Hs_n")
#check 3
test_that("mouse and human genes are different",{
expect_false(isTRUE(all.equal(bmt$Mm_n,bmt$Hs_n)))})
#continue
bmt$Mm_n <- gsub(bmt$Mm_n, pattern = "mir", replacement = "miR")
bmt$Hs_n <- gsub(bmt$Hs_n, pattern = "mir", replacement = "miR")
mh_m <- merge(x = miR, y = bmt, by.x = "Gene", by.y = "Mm_n", all = TRUE) 
mh_d <- mh_m[!duplicated(mh_m$name),] 
rownames(mh_d) <- mh_d$name
mh_d$Hs_n[is.na(mh_d$Hs_n)] <- as.character(mh_d$Gene[is.na(mh_d$Hs_n)])
mh_d$Hs_n <- gsub(mh_d$Hs_n, pattern = 'mmu', replacement = 'hsa')
#check 4
test_that("aspects of nodups are as expected", {
expect_equal(length(rownames(mh_d)),
length(rownames(miR)))
expect_equal(length(names(mh_d)),
length(names(mh_m)))
})
#continue
mh_d[order(mh_d$name),] -> mh_o
MicroRNA_full(mh_o$Hs_n, 'hsa') -> mh_o$microRNA
bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENTREZID",
OrgDb = org.Hs.eg.db) -> Y
Y[! duplicated(Y$GENENAME),] -> Y2
merge(x =mh_o, y=Y2, by.x="microRNA", by.y="GENENAME", all=TRUE)->miR_en
bitr(geneID = mh_o$microRNA, fromType = "GENENAME", toType = "ENSEMBL",
OrgDb = org.Hs.eg.db) -> Z
Z[! duplicated(Z$GENENAME),] -> Z2
merge(x=miR_en, y=Z2, by.x="microRNA", by.y="GENENAME", all=TRUE) -> miR_IDs
#check 5
test_that("miR_IDs should have entrezID and Ensembl data", {
expect_equal(length(names(miR_IDs)), 16)
expect_equal(length(rownames(miR_IDs)), length(rownames(miR)))
})
#continue
non_unique(Col = miR_IDs$Hs_n, sep = '-', suffix = 'p') -> miR_IDs$Hs_n
gsub(miR_IDs$Hs_n, pattern = "-1p", replacement = "-3p") -> miR_IDs$Hs_n
gsub(miR_IDs$Hs_n, pattern = "-2p", replacement = "-5p") -> miR_IDs$Hs_n
miR_IDs$ENTREZID -> miR_IDs$ENTREZID_adj
non_unique(Col=miR_IDs$ENTREZID_adj,sep ='.',suffix='')->miR_IDs$ENTREZID_adj
miR_IDs$ENSEMBL -> miR_IDs$ENSEMBL_adj
non_unique(Col=miR_IDs$ENSEMBL_adj,sep ='.',suffix='')->miR_IDs$ENSEMBL_adj
#check 6
test_that("new columns in miR_IDs are different from input columns", {
expect_false(isTRUE(all.equal(miR_IDs$ENTREZID, miR_IDs$ENTREZID_adj)))
expect_false(isTRUE(all.equal(miR_IDs$ENSEMBL, miR_IDs$ENSEMBL_adj)))
expect_equal(length(rownames(miR_IDs)), length(rownames(miR)))
})
#continue
MAE2 <- MultiAssayExperiment()
MAE2@ExperimentList$miR_entrez <- as.data.frame(cbind(GENENAME=rownames(
miR_IDs),ID=miR_IDs$ENTREZID))
MAE2@ExperimentList$miR_ensembl <- as.data.frame(cbind(GENENAME=rownames(
miR_IDs),ID=miR_IDs$ENSEMBL))
MAE2@ExperimentList$miR_adjusted_entrez <- as.data.frame(cbind(
GENENAME = rownames(miR_IDs), ID = miR_IDs$ENTREZID_adj))
MAE2@ExperimentList$miR_adjusted_ensembl <- as.data.frame(cbind(
GENENAME = rownames(miR_IDs), ID = miR_IDs$ENSEMBL_adj))
miR_IDs$microRNA<-miR_IDs$Gene<-miR_IDs$name<-miR_IDs$Hs_n<-miR_IDs$ENTREZID<-
miR_IDs$ENSEMBL<-miR_IDs$ENTREZID_adj<-miR_IDs$ENSEMBL_adj<-NULL
MAE2@ExperimentList$miR_human_renamed <- miR_IDs
#check 7
test_that("manual and functional outputs should be the same", {
expect_equal(dim(MAE@ExperimentList$miR_ensembl),
dim(MAE2@ExperimentList$miR_ensembl))
expect_equal(dim(MAE@ExperimentList$miR_adjusted_ensembl), 
dim(MAE2@ExperimentList$miR_adjusted_ensembl))
expect_equal(dim(MAE@ExperimentList$miR_entrez), 
dim(MAE2@ExperimentList$miR_entrez))
expect_equal(dim(MAE@ExperimentList$miR_adjusted_entrez),
dim(MAE2@ExperimentList$miR_adjusted_entrez))
expect_equal(dim(MAE@ExperimentList$miR_human_renamed), 
dim(MAE2@ExperimentList$miR_human_renamed))
})

