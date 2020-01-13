setwd("~/Documents/Package/smiRk/tests/testthat/")
#devtools::uses_testthat()
library(smiRk)
library(testthat)
#load filtered_genelist
readRDS("AddIDs_s/ensembl_genes.rds") -> ensembl_genes
#test function
eNames(method = 's', gene_IDs = ensembl_genes, ID_Column = 4) -> e_list
#check 1
test_that("e_list has 10 lists", {
expect_equal(class(e_list), "list")
expect_equal(length(e_list), 10)
})
#internal checks
sapply(ensembl_genes, function(x){sapply(x, `[[`, 4)}) -> Y
#check 2
test_that("Y is a matrix of 10 lists long", {
expect_equal(class(Y), "matrix")
expect_equal(length(Y), 10)
})
#continue
vapply(ensembl_genes, function(x){list(names(x))}, FUN.VALUE = list(1)) -> X
#check 3
test_that("X is a list of 2 lists", {
expect_equal(class(X), "list")
expect_equal(length(X), 2)
expect_equal(length(X[[1]]), 5)
})
#continue
unlist(X) -> Xnames
names(Y) <- Xnames
lapply(Y, function(x){ x[complete.cases(x)]}) -> y
#check 4
test_that("outputs are the same", {
expect_identical(e_list, y)
})
#save data
saveRDS(e_list, file = "eNames_s/elist_s.rds", compress = "xz")
