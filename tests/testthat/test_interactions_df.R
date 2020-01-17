#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load data
readRDS("miR_log2fc.rds") -> miR_log2fc
miR_log2fc[1:100,] -> miR_log2fc
readRDS("mRNA_log2fc.rds") -> mRNA_log2fc
mRNA_log2fc[1:200,] -> mRNA_log2fc
readRDS("test_net.rds") -> wiki_interest
#find which mRNAs are found in the wikipathway of interest
mRNA_log2fc[which(mRNA_log2fc$ID %in% wiki_interest$gene == TRUE),] -> mRNA_wiki
#check 1
#should be less in mRNA_wiki than mRNA_log2fc
test_that("mRNA_wiki is smaller than mRNA_log2fc", {
expect_lt(length(rownames(mRNA_wiki)),
length(rownames(mRNA_log2fc)))
})
name_combn <- expand.grid(rownames(miR_log2fc), rownames(mRNA_wiki),
stringsAsFactors = FALSE)
interactions_df <- do.call(rbind,apply(name_combn, 1, CorrTable,
miR_exprs = miR_log2fc,
mRNA_wiki = mRNA_wiki, maxInt = 5))
#continue
#check 2
#aspects of interactions_df
test_that("aspects of interactions_df are as expected", {
expect_equal(as.integer(length(rownames(miR_log2fc))*
length(rownames(mRNA_wiki))),
as.integer(length(rownames(interactions_df))))
expect_equal(length(names(interactions_df)),
5)
})
#save data
saveRDS(interactions_df, "interactions.rds",
compress = "xz")

