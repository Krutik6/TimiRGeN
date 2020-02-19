#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(org.Mm.eg.db)
library(clusterProfiler)
library(MultiAssayExperiment)
#load data
readRDS("IDs_mouse_miR.rds") -> miR
readRDS("filt_df.rds") -> Filt_df
MakeMapp(filt_df = Filt_df, 
miR_IDs_adj = miR@ExperimentList$miR_adjusted_entrez,
Datatype = 'L') -> MAPPdata
#check 1
#check aspects of MAPPdata
test_that("MAPPdata runs as expected", {
expect_equal(length(rownames(Filt_df)),
length(rownames(MAPPdata)))
expect_equal(length(names(MAPPdata)), 3)
})
#check 2
#check if 'En' works
MakeMapp(filt_df = Filt_df, 
miR_IDs_adj = miR@ExperimentList$miR_adjusted_ensembl,
Datatype = 'En') -> MAPPdata2
test_that("En and L lead to different results", {
expect_equal(MAPPdata[,1], MAPPdata2[,1])
expect_false(isTRUE(all.equal(MAPPdata, MAPPdata2)))
})


