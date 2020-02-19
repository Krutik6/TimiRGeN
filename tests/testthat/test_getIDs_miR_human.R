#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)
library(MultiAssayExperiment)
#load data
hs_miR -> miR
miR[1:20,] -> miR
#keep in vignette to show that the user may have to alter their input files
gsub(rownames(miR), pattern = "\\.", replacement =  "-") -> rownames(miR)
StartObject(miR = miR, mRNA = NULL) -> MAE
#test function
getIDs_miR_human(MAE, MAE@ExperimentList$miR) -> MAE
MAE@ExperimentList$miR_ensembl -> miR_ensembl
MAE@ExperimentList$miR_adjusted_ensembl -> miR_adjusted_ensembl
MAE@ExperimentList$miR_entrez -> miR_entrez
MAE@ExperimentList$miR_adjusted_entrez -> miR_adjusted_entrez
#check 1
test_that("output files should have two columns and the same number of
rownames", {
expect_equal(length(names(miR_ensembl)), 2)
expect_equal(length(names(miR_adjusted_ensembl)), 2)
expect_equal(length(names(miR_entrez)), 2)
expect_equal(length(names(miR_adjusted_entrez)), 2)
expect_equal(names(miR_adjusted_ensembl), names(miR_adjusted_entrez))
expect_equal(length(rownames(miR_adjusted_entrez)),
length(rownames(miR_ensembl)))
})
#internal checks
miR$Genes <- miR$MicroRNA <- rownames(miR)
gsub(x = miR$MicroRNA, pattern = "-3p", replacement = "") -> miR$MicroRNA
gsub(x = miR$MicroRNA, pattern = "-5p", replacement = "") -> miR$MicroRNA
#check 2
test_that("miR has 8 columns and the Genes and MicroRNA columns are different",
{
expect_equal(length(names(miR)), 8)
})

#continue
MicroRNA_full(miRdf = miR$MicroRNA, species = 'hsa') -> miR$MicroRNA
bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', toType = 'ENTREZID',
OrgDb = org.Hs.eg.db) -> miR_entrez_manual
bitr(geneID = miR$MicroRNA, fromType = 'GENENAME', toType = 'ENSEMBL',
OrgDb = org.Hs.eg.db) -> miR_ensembl_manual
merge(x = miR, y = miR_ensembl_manual, by.x = 'MicroRNA', by.y = 'GENENAME',
all = TRUE) -> miR_merged
merge(x = miR_merged, y = miR_entrez_manual, by.x = 'MicroRNA',
by.y = 'GENENAME', all = TRUE) -> miR_merged
miR_merged[!duplicated(miR_merged$Genes),] -> miR_merged
miR_merged[order(miR_merged$Genes),] -> miR_merged
#check 3
test_that("miR_merged has similiar entries to miR", {
expect_equal(length(rownames(miR_merged)),
length(rownames(miR)))
miR_merged <- miR_merged[order(miR_merged$Genes),]
expect_equal( miR_merged$Genes, rownames(miR))
})
#continue
non_unique(Col = miR_merged$ENTREZID, sep = ".",
suffix = "") -> miR_merged$ENTREZID_adjusted
non_unique(Col = miR_merged$ENSEMBL, sep = ".",
suffix = "") -> miR_merged$ENSEMBL_adjusted
miR_merged[! duplicated(miR_merged$Genes),] -> miR_merged
miR_merged <- miR_merged[order(miR_merged$Genes),]
rownames(miR_merged) <- miR_merged$Genes
#check 4
test_that("miR_merged has adjusted information", {
expect_equal(length(names(miR_merged)), 12)
expect_equal(names(miR_merged)[11], "ENTREZID_adjusted")
expect_equal(names(miR_merged)[12], "ENSEMBL_adjusted")
})
#continue
MAE2 <- MultiAssayExperiment()
MAE2@ExperimentList$miR_entrez <- as.data.frame(cbind(GENENAME = 
rownames(miR_merged), ID = miR_merged$ENTREZID))
MAE2@ExperimentList$miR_adjusted_entrez <- as.data.frame(cbind(
GENENAME = rownames(miR_merged),ID = miR_merged$ENTREZID_adjusted))
MAE2@ExperimentList$miR_ensembl <- as.data.frame(cbind( GENENAME = 
rownames(miR_merged), ID = miR_merged$ENSEMBL))
MAE2@ExperimentList$miR_adjusted_ensembl <- as.data.frame(cbind( 
GENENAME = rownames(miR_merged), ID = miR_merged$ENSEMBL_adjusted))
#check 5
test_that("manual and functional output are the same", {
expect_equal(MAE2@ExperimentList$miR_ensembl, 
MAE@ExperimentList$miR_ensembl)
expect_equal(MAE2@ExperimentList$miR_adjusted_ensembl,
MAE@ExperimentList$miR_adjusted_ensembl)
expect_equal(MAE2@ExperimentList$miR_entrez,
MAE@ExperimentList$miR_entrez)
expect_equal(MAE2@ExperimentList$miR_adjusted_entrez,
MAE@ExperimentList$miR_adjusted_entrez)
})
#save output
saveRDS(MAE, "IDs_human_miR.rds", compress = "xz")
