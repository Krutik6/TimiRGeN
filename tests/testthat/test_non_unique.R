#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load filtered_genelist
miR_IDs <- readRDS("IDs_mouse_miR.rds")
X <- assay(miR_IDs, 3)

miR_IDs_dup <- non_unique(Col = X$ID, sep = '-', suffix = 'p')

expect_false(isTRUE(all.equal(X$ID,  miR_IDs_dup)))
