#devtools::uses_testthat()
library(TimiRGeN)
library(testthat)
library(org.Mm.eg.db)
#load wpid2pathway data
MAE <- readRDS("wpdata.rds")
path_gene <- MAE[[1]]
path_data <- MAE[[3]]
#test GMT_ensembl function
MAE <- gmtEnsembl(MAE = MAE, path_gene = path_gene,
                   path_data = path_data, orgDB = org.Mm.eg.db)
#internal checks
ensembl_values <- suppressWarnings(bitr(geneID = path_gene$gene,
                                        fromType = 'ENTREZID',
                                        toType = 'ENSEMBL',
                                        OrgDb = org.Mm.eg.db))
#check 1
#ensembl values and entrezIDs
test_that("ensembl_values has 2 columns", {
    expect_that(length(names(ensembl_values)), equals(2))
})
#continue
path_gene_ensembl_manual <- merge(x = path_gene, y = ensembl_values,
                                  by.x = 'gene', by.y = 'ENTREZID')
#check 2
#wpid2gene_ensembl_manual have three columns
test_that("wpid2gene_ensembl_manual has 2 columns", {
    expect_that(length(names(path_gene_ensembl_manual)), equals(3))
})
#continue
path_data_ensembl_manual <- merge(x = path_gene_ensembl_manual,
                                  y = path_data, by.x = 'gene',
                                  by.y = 'gene')
#check 3
#wp_data_ensembl_manual have five columns
test_that("wp_data_ensembl_manual has 5 columns", {
    expect_that(length(names(path_data_ensembl_manual)), equals(5))
})
#continue
path_data_ensembl_corrected <- path_data_ensembl_manual[
    path_data_ensembl_manual[,2] == path_data_ensembl_manual[,4],]
#check 4
#column 2 and 4 are the same
test_that("Columns 2 and 4 of DFs are the same but were different", {
    expect_equal(path_data_ensembl_corrected[,2],
                 path_data_ensembl_corrected[,4])
    expect_false(isTRUE(path_data_ensembl_manual[,2]),
                 equals(path_data_ensembl_manual[,4]))
})
#continue
path_gene_ensembl_manual$gene <- NULL
path_data_ensembl_corrected$gene <- path_data_ensembl_corrected$wpid.y <- NULL
names(path_gene_ensembl_manual)[2] <- 'gene'
names(path_data_ensembl_corrected)[2] <- 'gene'
#check 5
#colnames are the same between dataframes

assay(MAE, 4)

test_that("Column features are as expected", {
    expect_equal(names(path_data_ensembl_corrected), names(assay(MAE, 5)))
    expect_equal(names(path_gene_ensembl_manual), names(assay(MAE, 4)))
    expect_equal(length(names(path_data_ensembl_corrected)), 3)
    expect_equal(length(names(path_gene_ensembl_manual)), 2)
})
#check 6
#There are no duplicates in wp_data_ensembl_manual
X <- path_data_ensembl_corrected[apply(path_data_ensembl_corrected, 1,
                                       function(x) length(unique(x))) > 1,]
test_that("X and wp_data_ensembl_manual are the same", {
    expect_equal(X, path_data_ensembl_corrected)
})
