#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
#load filtered_genelist
readRDS("ensembl_genes_s.rds") -> ensembl_genes
#test function
eNames(method = 's', gene_IDs = ensembl_genes, ID_Column = 4) -> e_list
#check 1
test_that("e_list has 10 lists", {
expect_equal(class(e_list), "list")
expect_equal(length(e_list), 10)
})
#internal checks
sapply(ensembl_genes, function(x){sapply(x, `[[`, 4)}) -> Y
#continue
vapply(ensembl_genes, function(x){list(names(x))}, FUN.VALUE = list(1)) -> X
unlist(X) -> Xnames
names(Y) <- Xnames
lapply(Y, function(x){ x[complete.cases(x)]}) -> y

