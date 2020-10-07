#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load data
Log2FC <- readRDS("log2fc.rds")

Ints <- readRDS("interactions.rds")

#find which mRNAs are found in the wikipathway of interest
mRNA_wiki <- assay(Log2FC, 3)[which(
             assay(Log2FC, 3)$ID %in% Ints$gene == TRUE),
             ]

#check 1
#should be less in mRNA_wiki than mRNA_log2fc
test_that("mRNA_wiki is smaller than mRNA_log2fc", {
    expect_lt(length(rownames(mRNA_wiki)),
    length(rownames(assay(Log2FC, 3))))
})

Log2FC

name_combn <- expand.grid(rownames(assay(Log2FC, 4)),
                          rownames(mRNA_wiki),
                          stringsAsFactors = FALSE)

interactions_df <- do.call(rbind,apply(name_combn, 1, corrTable,
                                       miR_exprs = assay(Log2FC, 4),
                                       mRNA_wiki = mRNA_wiki, maxInt = 5))

#continue
#check 2
#aspects of interactions_df
test_that("aspects of interactions_df are as expected", {
    expect_equal(length(names(interactions_df)), 5)
})

#save data
saveRDS(interactions_df, "corrmat.rds",compress = "xz")
