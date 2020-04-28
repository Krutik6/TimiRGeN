#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
#load filtered_genelist
ensembl_genes <- readRDS("ensembl_genes_s.rds")
MAE <- MultiAssayExperiment()
metadata(MAE)[["ensembl_genes"]] <- ensembl_genes
#test function
MAE<- eNames(MAE, method = 's',  gene_IDs = metadata(MAE)[[1]],
              ID_Column = 4)
#check 1
test_that("e_list has 10 lists", {
    expect_equal(class(metadata(MAE)[[2]]), "list")
    expect_equal(length(metadata(MAE)[[2]]), 10)
})
#internal checks
Y <- sapply(ensembl_genes, function(x){sapply(x, `[[`, 4)})
#continue
X <- vapply(ensembl_genes, function(x){list(names(x))}, FUN.VALUE = list(1))
Xnames <- unlist(X)
names(Y) <- Xnames
ID_list <- lapply(Y, function(x){ x[complete.cases(x)]})
metadata(MAE)[["ID_list"]] <- ID_list
