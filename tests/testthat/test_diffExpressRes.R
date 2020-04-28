#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load data
#format miR and mRNA data
#load ID data
MAE <- readRDS("MAE_mm.rds")
miR <- readRDS("IDs_mouse_miR.rds")
mRNA <- readRDS("IDs_mouse_mRNA.rds")
#test function
MAE <- diffExpressRes(MAE, df = assay(MAE, 2), dataType = 'Log2FC',
                       genes_ID = assay(mRNA, 1),
                       idColumn = 'GENENAME',
                       'mRNA_log2fc')
#internal checks on mRNA data
exp <- cbind(names = rownames(assay(MAE, 2)),
             assay(MAE, 2)[,grep('Log2FC', colnames(assay(MAE, 2)))])
#check 1
#exp should be the same row length as mRNA
test_that("aspects of exp should be similiar to mRNA", {
    expect_equal(length(rownames(assay(MAE, 2))), length(rownames(exp)))
    expect_equal(rownames(assay(MAE, 2)), rownames(exp))
})
#continue
merged <- merge(exp, mRNA@ExperimentList$mRNA_entrez, by.x = 'names',
                by.y = 'GENENAME',all = TRUE)
rownames(merged) <- merged[[1]]
merged[[1]] <- NULL
#check 2
#functional and manual output should be the same
test_that("output is the same in merged and mRNA_log2fc", {
    expect_equal(merged, assay(MAE, 3))
})
#Express on miRNA data
MAE <- diffExpressRes(MAE, df = assay(MAE, 1),
                    dataType = 'Log2FC',
                    genes_ID = assay(miR, 3),
                    idColumn = 'GENENAME',
                    name = 'miR_log2fc')
#save data as rds
saveRDS(MAE, "log2fc.rds", compress = "xz")
