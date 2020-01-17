#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
#format miR and mRNA data
mm_mRNA -> mRNA
mRNA[1:200,] -> mRNA
mm_miR -> miR
miR[1:100,] -> miR
#load ID data
readRDS("mm_miR_entrez.rds") -> miR_entrez
readRDS("mm_mRNA_entrez.rds") -> mRNA_entrez
#test function
Express(df = mRNA, dataType = 'Log2FC', genes_ID = mRNA_entrez,
idColumn = 'GENENAME') -> mRNA_log2fc
#internal checks on mRNA data
exp <- cbind(names = rownames(mRNA), mRNA[,grep('Log2FC', colnames(mRNA))])
#check 1
#exp should be the same row length as mRNA
test_that("aspects of exp should be similiar to mRNA", {
expect_equal(length(rownames(mRNA)),
length(rownames(exp)))
expect_equal(rownames(mRNA), rownames(exp))
})
#continue
merge(exp, mRNA_entrez, by.x = 'names', by.y = 'GENENAME',
all = TRUE) -> merged
rownames(merged) <- merged[[1]]
merged[[1]] <- NULL
#check 2
#rownames should be the same in merged and mRNA
test_that("rownames in merged are the same as in mRNA", {
expect_equal(length(rownames(mRNA)),
length(rownames(merged)))
mRNA <- mRNA[order(rownames(mRNA)),]
expect_equal(rownames(mRNA), rownames(merged))
})
#check 3
#functional and manual output should be the same
test_that("output is the same in merged and mRNA_log2fc", {
expect_equal(merged, mRNA_log2fc)
})

#Express on miRNA data
Express(df = miR, dataType = 'Log2FC', genes_ID = miR_entrez,
idColumn = 'GENENAME') -> miR_log2fc
#save data as rds
saveRDS(miR_log2fc, "miR_log2fc.rds", compress = "xz")
saveRDS(mRNA_log2fc, "mRNA_log2fc.rds", compress = "xz")
