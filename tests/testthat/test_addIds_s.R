#devtools::uses_testthat()
library(TimiRGeN)
library(MultiAssayExperiment)
library(testthat)
#load filtered_genelist
MAE <- MultiAssayExperiment()
metadata(MAE)[["filtered_genelist"]] <- readRDS("filtered_genelist_s.rds")
#load mouse data
miR <- readRDS("IDs_mouse_miR.rds")
mRNA <- readRDS("IDs_mouse_mRNA.rds")
#test function
MAE <- addIds(MAE, method = 's',
            filtered_genelist = metadata(MAE)[[1]],
            miR_IDs = assay(miR, 2),
            mRNA_IDs = assay(mRNA, 2))
#check 1
test_that("ensembl_genes has qualities which are expected", {
    expect_equal(class(metadata(MAE)[[2]]), "list")
    expect_equal(length(metadata(MAE)[[1]]), length(metadata(MAE)[[2]]))
})
#continue
metadata(MAE)[["data_IDs2"]] <- Map(function(x, y) lapply(x, function(dat) {
                                     dat$GENENAME <- row.names(dat);
                                     merge(dat, y)}), metadata(MAE)[[1]],
                                     list(assay(miR, 2), assay(mRNA, 2)))
#check 2
test_that("outputs are the same", {
    expect_equal(metadata(MAE)[[3]], metadata(MAE)[[2]])
})

#save data
saveRDS(metadata(MAE)[[2]],file = "ensembl_genes_s.rds", compress = "xz")

