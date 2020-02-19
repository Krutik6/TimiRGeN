#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
#load filtered_genelist
MAE <- MultiAssayExperiment()
MAE@metadata$filtered_genelist <- readRDS(
"filtered_genelist_s.rds") 
#load mouse data
miR <- readRDS("IDs_mouse_miR.rds")
mRNA <- readRDS("IDs_mouse_mRNA.rds")
#test function
AddIDs(method = 's', filtered_genelist = MAE@metadata$filtered_genelist, 
miR_IDs = miR@ExperimentList$miR_ensembl, 
mRNA_IDs = mRNA@ExperimentList$mRNA_ensembl) -> MAE@metadata$ensembl_genes
#check 1
test_that("ensembl_genes has qualities which are expected", {
expect_equal(class(MAE@metadata$ensembl_genes), "list")
expect_equal(length(MAE@metadata$ensembl_genes[[1]]),
length(MAE@metadata$filtered_genelist[[1]]))
expect_gt(length(MAE@metadata$ensembl_genes[[1]][[1]]),
length(MAE@metadata$filtered_genelist[[1]][[1]]))
})
#internal checks
colnames(miR@ExperimentList$miR_ensembl) <- c("GENENAME", "ID")
colnames(mRNA@ExperimentList$mRNA_ensembl) <- c("GENENAME", "ID")

Map(function(x, y) lapply(x, function(dat) {dat$GENENAME <- row.names(dat);
merge(dat, y)}), MAE@metadata$filtered_genelist,
list(miR@ExperimentList$miR_ensembl, mRNA@ExperimentList$mRNA_ensembl)) -> X
#check 2
test_that("outputs are the same", {
expect_equal(X, MAE@metadata$ensembl_genes)
})
#save data
saveRDS(MAE@metadata$ensembl_genes,
file = "ensembl_genes_s.rds", compress = "xz")
