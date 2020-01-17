#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load mouse data
readRDS("mm_miR_entrez.rds") -> miR_entrez
readRDS("mm_mRNA_entrez.rds") -> mRNA_entrez
#load gene_entrez
readRDS("gene_entrez.rds") -> gene_entrez
#perform function
eNames(method = 'c', gene_IDs = gene_entrez, ID_Column = 4) -> e_list
#internal checks
lapply(gene_entrez, function(x){x[[4]]}) -> e
lapply(e, function(x){ x[complete.cases(x)]}) -> y
#check 1
#manual output is the same as funcitonal output
test_that("Output is as expects", {
expect_equal(y, e_list)
})
#save data
saveRDS(e_list, "elist.rds", compress = "xz")
