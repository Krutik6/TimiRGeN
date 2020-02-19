#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load data
#format miR and mRNA data
#load ID data
readRDS("MAE_mm.rds") -> MAE
readRDS("IDs_mouse_miR.rds") -> miR
readRDS("IDs_mouse_mRNA.rds") -> mRNA
#test function
Express(df = MAE@ExperimentList$mRNA, dataType = 'Log2FC',
genes_ID = mRNA@ExperimentList$mRNA_entrez,
idColumn = 'GENENAME') -> MAE@ExperimentList$mRNA_log2fc
#internal checks on mRNA data
exp <- cbind(names = rownames(MAE@ExperimentList$mRNA),
MAE@ExperimentList$mRNA[,grep('Log2FC', colnames(MAE@ExperimentList$mRNA))])
#check 1
#exp should be the same row length as mRNA
test_that("aspects of exp should be similiar to mRNA", {
expect_equal(length(rownames(MAE@ExperimentList$mRNA)),
length(rownames(exp)))
expect_equal(rownames(MAE@ExperimentList$mRNA), rownames(exp))
})
#continue
merge(exp, mRNA@ExperimentList$mRNA_entrez, by.x = 'names', by.y = 'GENENAME',
all = TRUE) -> merged
rownames(merged) <- merged[[1]]
merged[[1]] <- NULL
#check 2
#functional and manual output should be the same
test_that("output is the same in merged and mRNA_log2fc", {
expect_equal(merged, MAE@ExperimentList$mRNA_log2fc)
})
#Express on miRNA data
Express(df = MAE@ExperimentList$miR, dataType = 'Log2FC', 
genes_ID = miR@ExperimentList$miR_entrez,
idColumn = 'GENENAME') -> MAE@ExperimentList$miR_log2fc
#save data as rds
saveRDS(MAE, "log2fc.rds", compress = "xz")
