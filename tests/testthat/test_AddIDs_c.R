#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load filtered_genelist
MAE <- MultiAssayExperiment()
MAE@metadata$filtered_genelist <- readRDS("filtered_genelist_c.rds") 
#load mouse data
miR <- readRDS("IDs_mouse_miR.rds")
mRNA <- readRDS("IDs_mouse_mRNA.rds")
#use function
AddIDs(method = 'c', filtered_genelist = MAE@metadata$filtered_genelist,
miR_IDs = miR@ExperimentList$miR_entrez, 
mRNA_IDs = mRNA@ExperimentList$mRNA_entrez) -> MAE@metadata$gene_entrez
#internal checks
colnames(miR@ExperimentList$miR_entrez) <- c("GENENAME", "ID")
colnames(mRNA@ExperimentList$mRNA_entrez) <- c("GENENAME", "ID")
#check 1
#colnames should be the same and the first column is GENENAME, second column
#is ID
test_that("miR_entrez and mRNA_entrez have the same column names",{
expect_equal(colnames(miR@ExperimentList$miR_entrez), 
colnames(mRNA@ExperimentList$mRNA_entrez))
expect_equal(colnames(miR@ExperimentList$miR_entrez[1]), "GENENAME")
expect_equal(colnames(miR@ExperimentList$miR_entrez[2]), "ID")
expect_equal(colnames(mRNA@ExperimentList$mRNA_entrez[1]), "GENENAME")
expect_equal(colnames(mRNA@ExperimentList$mRNA_entrez[2]), "ID")
})
#continue
rbind(miR@ExperimentList$miR_entrez,
mRNA@ExperimentList$mRNA_entrez) -> geneIDs
geneIDs[! duplicated(geneIDs$GENENAME),] -> genes_id
lapply(MAE@metadata$filtered_genelist, function(x){cbind(
'GENENAME' = rownames(x),
x)})-> X
#check 2
#each DF should have 3 columns
test_that("Each nested dataframe in X should have 1 more column
than input", {
for (i in 1:5) {
expect_gt(length(colnames(X[[i]])),
length(colnames(MAE@metadata$filtered_genelist[[i]])))
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
expect_equal(MAE@metadata$gene_entrez, Y)
})
#save data
saveRDS(MAE@metadata$gene_entrez, "gene_entrez_c.rds",
compress = "xz")

