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
df2 <- assay(MAE, 4) %>% mutate(TargetScan = as.integer(
                        rownames(assay(MAE, 4)) %in% assay(MAE, 1)[[1]]))
rownames(df2) <- rownames(assay(MAE, 4))
df3 <-df2 %>% mutate(miRDB = as.integer(rownames(df2) %in% assay(MAE, 2)[[1]]))
rownames(df3) <- rownames(assay(MAE, 4))
df4 <- df3 %>% mutate(Predicted_Interactions = TargetScan + miRDB)
rownames(df4) <- rownames(assay(MAE, 4))
df5 <- df4 %>% mutate(miRTarBase = as.integer(
               rownames(df4) %in% assay(MAE, 3)[[1]]))
rownames(df5) <- rownames(assay(MAE, 4))
df6 <- df5 %>% mutate(Pred_Fun = Predicted_Interactions + miRTarBase)
rownames(df6) <- rownames(assay(MAE, 4))
#check 2
#dfs should have increasing numbers of columns
test_that("there are increasing numbers of columns", {
    expect_gt(length(names(df2)), length(names(assay(MAE, 4))))
    expect_gt(length(names(df3)), length(names(df2)))
    expect_gt(length(names(df4)), length(names(df3)))
    expect_gt(length(names(df5)), length(names(df4)))
    expect_gt(length(names(df6)), length(names(df5)))
    expect_equal(length(names(df6)),10)
})
#check 3
#outcomes should be the same
test_that("functional and manual output is the same", {
    expect_equal(rownames(df6), rownames(assay(MAE, 5)))
    expect_equal(colnames(df6), colnames(assay(MAE, 5)))
})
#save data
saveRDS(df6, "MiningMatrix.rds", compress = "xz")
