#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load filtered_genelist
readRDS("IDs_human_miR.rds") -> miR_IDs

non_unique(Col = miR_IDs@ExperimentList$miR_entrez$ID, 
sep = '-', suffix = 'p') -> miR_IDs_dup

expect_false(isTRUE(all.equal(miR_IDs@ExperimentList$miR_entrez$ID, 
miR_IDs_dup)))
