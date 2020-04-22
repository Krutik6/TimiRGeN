#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(MultiAssayExperiment)
#load filtered_genelist
MAE <- MultiAssayExperiment()
metadata(MAE)[["filtered_genelist"]] <- readRDS("filtered_genelist_c.rds")
#load mouse data
miR <- readRDS("IDs_mouse_miR.rds")
mRNA <- readRDS("IDs_mouse_mRNA.rds")
#use function
MAE <- addIds(MAE, method = 'c',
              filtered_genelist = metadata(MAE)[[1]],
              miR_IDs = assay(miR, 1),
              mRNA_IDs = assay(mRNA, 1))

#internal checks
#check 1
#colnames should be the same and the first column is GENENAME, second column
#is ID
test_that("miR_entrez and mRNA_entrez have the same column names",{
    expect_equal(colnames(assay(miR, 1)), colnames(assay(mRNA, 1)))
    expect_equal(colnames(assay(miR, 1)[1]), "GENENAME")
    expect_equal(colnames(assay(miR, 1)[2]), "ID")
    expect_equal(colnames(assay(mRNA, 1)[1]), "GENENAME")
    expect_equal(colnames(assay(mRNA, 1)[2]), "ID")
})
#continue
geneIds <- rbind(assay(miR, 1), assay(mRNA, 1))
genes_id <- geneIds[! duplicated(geneIds$GENENAME),]
X <- lapply(metadata(MAE)[[1]], function(x){cbind('GENENAME' = rownames(x),x)})
#check 2
#each DF should have 3 columns
test_that("Each nested dataframe in X should have 1 more columnthan input", {
    for (i in 1:5) {
        expect_gt(length(colnames(X[[i]])),
        length(colnames(metadata(MAE)[[1]][[i]])))
    }
    })
#continue
Y <- lapply(X, function(x){merge(x, genes_id)})
#check 3
#check Y columns and column names
test_that("Y columns and column names", {
    for (i in 1:5) {
        expect_gt(length(colnames(Y[[i]])),
        length(colnames(X[[i]])))
        expect_true(colnames(Y[[i]][1]) == "GENENAME")
        expect_true(colnames(Y[[i]][4]) == "ID")
    }
})
#check 4
#check functional and manual outputs are the same
test_that("Functional output and manual outputs are the same",{
    expect_equal(metadata(MAE)[[2]], Y)
})
metadata(MAE)[["data_IDs"]] <- Y
#save data
saveRDS(metadata(MAE)[[2]], "gene_entrez_c.rds", compress = "xz")
