#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
#load Data
Clusters <- readRDS("Clusters.rds")
#test function
Clusters <- ReturnCluster(MAE = Clusters, clusterData = assay(Clusters, 2), 
                          whichCluster = 1, fitCluster = 0.1)

#check 1
test_that("returned dataframe is smaller than original", {
    expect_lt(length(rownames(assay(Clusters, 4))),
    length(rownames(assay(Clusters, 2))))
})
#continue
Clusters <- ReturnCluster(MAE = Clusters, clusterData = assay(Clusters, 2), 
                          whichCluster = 1, fitCluster = 0.5)

a <- "Cluster:" 
cl <- 1
b <- "_fit:"
fit <- 0.1


