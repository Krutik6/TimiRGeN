setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(clusterProfiler)
library(org.Mm.eg.db)
# load data
mm_miR -> miR
miR[1:50,]-> miR
#test function
getIDs_miR_mouse(miR)
#test outputs have expected number of columns
test_that("miR_Id data have two columns", {
expect_that(as.numeric(ncol(miR_entrez)), equals(2))
expect_that(as.numeric(ncol(miR_ensembl)), equals(2))
expect_that(as.numeric(ncol(miR_adjusted_entrez)), equals(2))
expect_that(as.numeric(ncol(miR_adjusted_ensembl)), equals(2))
})
#internal checks
miR$Genes <- miR$names <- rownames(miR)
sub(x = miR$Genes, pattern = "-3p", replacement = "") -> miR$Genes
sub(x = miR$Genes, pattern = "-5p", replacement = "") -> miR$Genes
#check 1
#test that -3p and -5p have been removed
test_that("-3p and -5p should be absent", {
expect_that(length(grep(miR$Genes, pattern = "-5p")), equals(0))
expect_that(length(grep(miR$Genes, pattern = "-3p")), equals(0))
})
#continue
MicroRNA_full(miRdf = miR$Genes, species = 'mmu') -> miR$Genes
#check 2
test_that("Genes and names are different", {
expect_false(isTRUE(all.equal(miR$Genes, miR$names)))
})
#continue
bitr(geneID = miR$Genes, fromType = 'GENENAME', toType = 'ENTREZID',
OrgDb = org.Mm.eg.db) -> miR_entrez2
bitr(geneID = miR$Genes, fromType = 'GENENAME', toType = 'ENSEMBL',
OrgDb = org.Mm.eg.db) -> miR_ensembl2
#check 3
test_that("miR_entrez and miR_ensembl are different but similiar", {
expect_equal(length(names(miR_ensembl2)), 2)
expect_equal(length(names(miR_entrez2)), 2)
})
#conitnue
merge(x = miR, y = miR_ensembl2, by.x = 'Genes', by.y = 'GENENAME',
all = TRUE) -> miR_merged
merge(x = miR_merged, y = miR_entrez2, by.x = 'Genes', by.y = 'GENENAME',
all = TRUE) -> miR_merged
miR_merged[!duplicated(miR_merged$names),] -> miR_merged
miR_merged[order(miR_merged$names),] -> miR_merged
#check 4
test_that("miR_merged has 14 columns", {
expect_equal(length(names(miR_merged)), 14)
})
#continue
ave(as.character(miR_merged$ENTREZID), miR_merged$ENTREZID,
FUN = function(x) if (length(x)>1)
paste0(x[1], '.', seq_along(x), '')
else x) -> miR_merged$ENTREZID_adjusted
ave(as.character(miR_merged$ENSEMBL), miR_merged$ENSEMBL,
FUN = function(x) if (length(x)>1)
paste0(x[1], '.', seq_along(x), '')
else x) -> miR_merged$ENSEMBL_adjusted
#check 4
test_that("miR_merged has 16 columns", {
expect_equal(length(names(miR_merged)), 16)
})
#continue
miR_entrez_manual <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENTREZID))
miR_adjusted_entrez_manual <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENTREZID_adjusted))
miR_ensembl_manual <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENSEMBL))
miR_adjusted_ensembl_manual <- as.data.frame(cbind(GENENAME = miR_merged$names,
ID = miR_merged$ENSEMBL_adjusted))
#check 5
test_that("manual and functional output are the same", {
expect_equal(miR_ensembl, miR_ensembl_manual)
expect_equal(miR_adjusted_ensembl, miR_adjusted_ensembl_manual)
expect_equal(miR_entrez, miR_entrez_manual)
expect_equal(miR_adjusted_entrez, miR_adjusted_entrez_manual)
})
#save output
setwd("getIDs_miR_mouse/")
saveRDS(miR_adjusted_ensembl, "miR_adjusted_ensembl.rds", compress = "xz")
saveRDS(miR_adjusted_entrez, "miR_adjusted_entrez.rds", compress = "xz")
saveRDS(miR_ensembl, "miR_ensembl.rds", compress = "xz")
saveRDS(miR_entrez, "miR_entrez.rds", compress = "xz")


