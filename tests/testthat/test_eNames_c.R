#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load gene_entrez
ID_names <- readRDS("gene_entrez_c.rds")

MAE <- MultiAssayExperiment()

metadata(MAE)[["ID_names"]] <- ID_names

#perform function
MAE <- eNames(MAE, method = 'c', gene_IDs = metadata(MAE)[[1]],
              ID_Column = 4)

#internal checks
e <- lapply(ID_names, function(x){x[[4]]})

y <- lapply(e, function(x){ x[complete.cases(x)]})

#check 1
#manual output is the same as functional output
test_that("Output is as expects", {
    expect_equal(length(e_list), 5)
    expect_equal(length(y), 5)
})

metadata(MAE)[["ID_list"]] <- y
