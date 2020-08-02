#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(org.Mm.eg.db)

#load data
MAE <- MultiAssayExperiment(list(
    TargetScans_res = as.data.frame(readRDS("TargetScans_results.rds")),
    miRDB_res = as.data.frame(readRDS("miRDB_resuts.rds")),
    miRTarBase_res = as.data.frame(readRDS("miRTarBase_results.rds")),
    Interactions = readRDS("corrmat.rds")
))

#test DataMiningMatrix
MAE <- dataMiningMatrix(MAE,
                      corrTable = assay(MAE, 4),
                      targetscan = assay(MAE, 1),
                      mirdb = assay(MAE, 2),
                      mirtarbase = assay(MAE, 3))

#check 1
#X should have 10 columns and the last 5 should be numeric
test_that("last 5 columns are numerics", {
    expect_equal(length(names(assay(MAE, 5))), 10)
    for (i in 6:10) {
        expect_true(is.numeric(assay(MAE, 5)[,i]))
}
})
#internal checks
X <- assay(MAE, 4)

X$Pred_Fun <- X$miRTarBase <- X$Pred_only <- X$miRDB <-
             X$TargetScan <-numeric(nrow(X))

X$TargetScan <- as.integer(rownames(X) %in% assay(MAE, 1)[1])


X$miRDB <- as.integer(rownames(X) %in% assay(MAE, 2)[1])

X$Pred_only <- X$TargetScan + X$miRDB

X$miRTarBase <- as.integer(rownames(X) %in% assay(MAE, 3)[1])

X$Pred_Fun <- X$Pred_only + X$miRTarBase

#check 2
#outcomes should be the same
test_that("functional and manual output is the same", {
    expect_equal(rownames(X), rownames(assay(MAE, 5)))
    expect_equal(colnames(X), colnames(assay(MAE, 5)))
})

#save data
saveRDS(X, "MiningMatrix.rds", compress = "xz")
