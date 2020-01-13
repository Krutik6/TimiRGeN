setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load filtered_genelist
readRDS("SignificantVals_c/filtered_genelist.rds") -> filtered_genelist
#load mouse data
readRDS("getIDs_miR_mouse/miR_entrez.rds") -> miR_entrez
readRDS("getIDs_mRNA_mouse/mRNA_entrez.rds") -> mRNA_entrez
#use function
AddIDs(method = 'c', filtered_genelist = filtered_genelist,
miR_IDs = miR_entrez, mRNA_IDs = mRNA_entrez) -> gene_entrez
#internal checks
colnames(miR_entrez) <- c("GENENAME", "ID")
colnames(mRNA_entrez) <- c("GENENAME", "ID")
#check 1
#colnames should be the same and the first column is GENENAME, second column
#is ID
test_that("miR_entrez and mRNA_entrez have the same column names",{
expect_equal(colnames(miR_entrez), colnames(mRNA_entrez))
expect_equal(colnames(miR_entrez[1]), "GENENAME")
expect_equal(colnames(miR_entrez[2]), "ID")
expect_equal(colnames(mRNA_entrez[1]), "GENENAME")
expect_equal(colnames(mRNA_entrez[2]), "ID")
})
#continue
rbind(miR_entrez, mRNA_entrez) -> geneIDs
geneIDs[! duplicated(geneIDs$GENENAME),] -> genes_id
lapply(filtered_genelist, function(x){cbind('GENENAME' = rownames(x),
x)})-> X
#check 2
#each DF should have 3 columns
test_that("Each nested dataframe in X should have 1 more column
than input", {
for (i in 1:5) {
expect_gt(length(colnames(X[[i]])),
length(colnames(filtered_genelist[[i]])))
}
})
#continue
lapply(X, function(x){merge(x, genes_id)}) -> Y
#check 3
#check Y columns and column names
test_that("Y columns and column names", {
for (i in 1:5) {
expect_gt(length(colnames(Y[[i]])),
length(colnames(X[[i]])))
expect_true(colnames(Y[[i]][1]) == "GENENAME")
expect_true(colnames(Y[[i]][4]) == "ID")
}
})
#check 4
#check functional and manual outputs are the same
test_that("Functional output and manual outputs are the same",{
expect_equal(gene_entrez, Y)
})
#save data
saveRDS(gene_entrez, "AddIDs_c/gene_entrez.rds", compress = "xz")

