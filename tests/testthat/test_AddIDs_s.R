#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load filtered_genelist
readRDS(file = "filtered_genelist_s") -> filtered_genelist
readRDS(file = "mm_miR_ensembl.rds") -> miR_ensembl
readRDS(file = "mm_mRNA_ensembl.rds") -> mRNA_ensembl
#test function
AddIDs(method = 's', filtered_genelist = filtered_genelist,
miR_IDs = miR_ensembl, mRNA_IDs = mRNA_ensembl) -> ensembl_genes
#check 1
test_that("ensembl_genes has qualities which are expected", {
expect_equal(class(ensembl_genes), "list")
expect_equal(length(ensembl_genes[[1]]), length(filtered_genelist[[1]]))
expect_gt(length(ensembl_genes[[1]][[1]]),
length(filtered_genelist[[1]][[1]]))
})
#internal checks
colnames(miR_ensembl) <- c("GENENAME", "ID")
colnames(mRNA_ensembl) <- c("GENENAME", "ID")
Map(function(x, y) lapply(x, function(dat) {dat$GENENAME <- row.names(dat);
merge(dat, y)}), filtered_genelist, list(miR_ensembl, mRNA_ensembl)) -> X
#check 2
test_that("outputs are the same", {
expect_equal(X, ensembl_genes)
})
#save data
saveRDS(ensembl_genes, file = "ensembl_genes.rds", compress = "xz")
