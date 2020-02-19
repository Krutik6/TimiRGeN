#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
library(Mfuzz)
#load data
readRDS("Clusters.rds") -> Clusters
#check function
Quickfuzz(Mfuzzdata = Clusters@ExperimentList$MfuzzData,
Clusters = Clusters@metadata$Clusters, W = FALSE)
dev.off()
#there are 5 graphs showing indication 5 clusters
test_that("Clusters and Mfuzzdata have 5 clusters", {
expect_equal(max(unique(Clusters@metadata$Clusters$cluster)), 5)
expect_equal(length(rownames(Clusters@ExperimentList$MfuzzData@phenoData)), 5)
})

