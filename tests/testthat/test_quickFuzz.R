#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)

#load data
Clusters <- readRDS("Clusters.rds")

#check function
quickFuzz(Mfuzzdata = experiments(Clusters)[[3]],
          Clusters = metadata(Clusters)[[1]], W = FALSE)

dev.off()

#there are 5 graphs showing indication 5 clusters
test_that("Clusters and Mfuzzdata have 5 clusters", {
    expect_equal(max(unique(Clusters@metadata$Clusters$cluster)), 5)
    expect_equal(length(rownames(Clusters@ExperimentList$MfuzzData@phenoData)),
                 5)
})
