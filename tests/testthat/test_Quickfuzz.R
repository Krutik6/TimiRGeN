#devtools::uses_testthat()
library(smiRk)
library(testthat)
library(Mfuzz)
#load data
readRDS("Clusters.rds") -> Clusters
readRDS("Mfuzzdata.rds") -> Mfuzzdata
#check function
Quickfuzz(Mfuzzdata = Mfuzzdata, Clusters = Clusters, W = FALSE)
dev.off()
#there are 5 graphs showing indication 5 clusters
test_that("Clusters and Mfuzzdata have 5 clusters", {
expect_equal(max(unique(Clusters$cluster)), 5)
expect_equal(length(rownames(Mfuzzdata@phenoData@data)), 5)
})

